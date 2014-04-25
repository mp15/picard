/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.sam.markduplicates;

import picard.PicardException;
import picard.sam.CoordinateSortedPairInfoMap;

import java.io.*;
import java.util.*;

/**
 * Disk-based implementation of ReadEndsMarkDuplicatesMap.  A subdirectory of the system tmpdir is created to store
 * files, one for each reference sequence.  The reference sequence that is currently being queried (i.e. the
 * sequence for which remove() has been most recently called) is stored in RAM.  ReadEnds for all other sequences
 * are stored on disk.
 *
 * When put() is called for a sequence that is the current one in RAM, the ReadEnds object is merely put into the
 * in-memory map.  If put() is called for a sequence ID that is not the current RAM one, the ReadEnds object is
 * appended to the file for that sequence, creating the file if necessary.
 *
 * When remove() is called for a sequence that is the current one in RAM, remove() is called on the in-memory map.
 * If remove() is called for a sequence other than the current RAM sequence, then the current RAM sequence is written
 * to disk, the new sequence is read from disk into RAM map, and the file for the new sequence is deleted.
 *
 * If things work properly, and reads are processed in genomic order, records will be written for mates that are in
 * a later sequence.  When the mate is reached in the input SAM file, the file that was written will be deleted.
 * This should result in all temporary files being deleted by the time all the reads are processed.  The temp
 * directory is marked to be deleted on exit so everything should get cleaned up.
 * 
 * @author alecw@broadinstitute.org
 */
class DiskReadEndsMarkDuplicatesMap implements ReadEndsMarkDuplicatesMap {
    private final CoordinateSortedPairInfoMap<String, ReadEndsMarkDuplicates> pairInfoMap;
    DiskReadEndsMarkDuplicatesMap(int maxOpenFiles) {
        pairInfoMap = new CoordinateSortedPairInfoMap<String, ReadEndsMarkDuplicates>(maxOpenFiles, new Codec());
    }

    public ReadEndsMarkDuplicates remove(int mateSequenceIndex, String key) {
        return pairInfoMap.remove(mateSequenceIndex, key);
    }

    public void put(int mateSequenceIndex, String key, ReadEndsMarkDuplicates readEnds) {
        pairInfoMap.put(mateSequenceIndex, key, readEnds);
    }

    public int size() {
        return pairInfoMap.size();
    }

    public int sizeInRam() {
        return pairInfoMap.sizeInRam();
    }

    private static class Codec implements CoordinateSortedPairInfoMap.Codec<String, ReadEndsMarkDuplicates> {
        private final ReadEndsMarkDuplicatesCodec readEndsMarkDuplicatesCodec = new ReadEndsMarkDuplicatesCodec();

        public void setInputStream(final InputStream is) {
            readEndsMarkDuplicatesCodec.setInputStream(is);
        }

        public void setOutputStream(final OutputStream os) {
            readEndsMarkDuplicatesCodec.setOutputStream(os);
        }

        public Map.Entry<String, ReadEndsMarkDuplicates> decode() {
            try {
                final String key = readEndsMarkDuplicatesCodec.getInputStream().readUTF();
                final ReadEndsMarkDuplicates record = readEndsMarkDuplicatesCodec.decode();
                return new AbstractMap.SimpleEntry<java.lang.String,ReadEndsMarkDuplicates>(key, record);
            } catch (IOException e) {
                throw new PicardException("Error loading ReadEndsMarkDuplicatesMap from disk", e);
            }
        }

        public void encode(final String key, final ReadEndsMarkDuplicates readEnds) {
            try {
                readEndsMarkDuplicatesCodec.getOutputStream().writeUTF(key);
                readEndsMarkDuplicatesCodec.encode(readEnds);
            } catch (IOException e) {
                throw new PicardException("Error spilling ReadEndsMarkDuplicatesMap to disk.", e);
            }
        }
    }
    
}