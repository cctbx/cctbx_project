from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format


class FormatBruker(Format):
    """cctbx has no authoritative sources describing the Bruker format.
    Positive identification:  when listing the header out in 16-byte chunks,
    there are numerous chunks of the form
    KEYNAME:xxxxxxxx
    KEYNAME:xxxxxxxx
    ...in other words a 7-upper case character (or space) keyname followed by
    a colon, starting at the beginning of every fifth chunk.
    Here we will take a series of 12 of these in the first 1024 characters
    as a positive fingerprint of Bruker format.
    """

    @staticmethod
    def read_header_lines(image_path, max_bytes=512 * 50):
        """Attempt to read the whole Bruker header into a list of lines of 80 char
        length. Finish reading lines when either the amount of data determined by
        HDRBLKS is consumed, or if max_bytes is reached, or a line contains a
        non-ASCII character. HDRBLKS is interpreted as being an integer number
        of 512 byte blocks. HDRBLKS must be divisible by 5, so that the total
        number of bytes in the header is divisible by 80."""

        hdr_lines = []
        max_bytes = (max_bytes // 80) * 80
        with FormatBruker.open_file(image_path, "rb") as f:
            while True:

                # read a line and look for HDRBLKS
                line = f.read(80)
                if not line:
                    break
                if line.startswith("HDRBLKS"):
                    try:
                        val = int(line.split(":", 1)[1])
                        if val % 5 == 0:
                            max_bytes = min(val * 512, max_bytes)
                    except ValueError:
                        pass

                try:
                    line.decode("ascii")
                except UnicodeDecodeError:
                    break

                # looks like a valid line
                hdr_lines.append(line)

                if max_bytes is not None:
                    if f.tell() >= max_bytes:
                        break

        return hdr_lines

    @staticmethod
    def parse_header(header_lines):
        header_dic = {}

        for l in header_lines:
            k_v = l.split(":", 1)
            if len(k_v) == 1:
                continue
            k, v = [v.strip() for v in k_v]
            if k in header_dic:
                header_dic[k] = header_dic[k] + "\n" + v
            else:
                header_dic[k] = v

        return header_dic

    @staticmethod
    def understand(image_file):
        try:
            tag = FormatBruker.open_file(image_file, "rb").read(1024)
        except IOError:
            return False
        matches = []
        for x in xrange(0, 1024, 80):
            word = tag[x : x + 16]
            if word[0:7].isupper() and word[7] == ":":
                matches.append(word)

        return len(matches) >= 12

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        Format.__init__(self, image_file, **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        """Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size."""
        from iotbx.detectors.bruker import BrukerImage

        self.detectorbase = BrukerImage(self._image_file)
        # self.detectorbase.readHeader() #unnecessary for the Bruker specialization

    def _goniometer(self):

        if self.detectorbase.parameters["OSC_RANGE"] > 0:
            return self._goniometer_factory.single_axis()
        else:
            return self._goniometer_factory.single_axis_reverse()

    def _detector(self):
        """Return a model for a simple detector"""

        twotheta = self.detectorbase.parameters["TWOTHETA"]
        # At present, ignore non-zero two theta for the dxtbx model
        # XXX Return to this issue later.

        return self._detector_factory.simple(
            sensor="CCD",
            distance=self.detectorbase.parameters["DISTANCE"],
            beam_centre=(
                self.detectorbase.parameters["BEAM_CENTER_X"],
                self.detectorbase.parameters["BEAM_CENTER_Y"],
            ),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(
                self.detectorbase.parameters["PIXEL_SIZE"],
                self.detectorbase.parameters["PIXEL_SIZE"],
            ),
            image_size=(
                self.detectorbase.parameters["SIZE1"],
                self.detectorbase.parameters["SIZE2"],
            ),
            trusted_range=(0, self.detectorbase.parameters["CCD_IMAGE_SATURATION"]),
            mask=[],
        )  # a list of dead rectangles

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.parameters["WAVELENGTH"])

    def _scan(self):
        """Return the scan information for this image."""

        return self._scan_factory.single(
            filename=self._image_file,
            format="BrukerCCD",
            # It's not at all clear how to recover the exposure time from the header
            # or even whether it is recorded.
            # XXX Here it will simply be set to a default number.
            exposure_times=1,
            osc_start=self.detectorbase.parameters["OSC_START"],
            osc_width=abs(self.detectorbase.parameters["OSC_RANGE"]),
            epoch=None,
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatBruker.understand(arg))
