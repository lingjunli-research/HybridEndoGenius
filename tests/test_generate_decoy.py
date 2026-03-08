"""
Unit tests for workflow/generate_decoy.py

Tests cover:
- Output file creation
- Spectrum count preservation
- Peak count preservation per spectrum
- Metadata/header line copying
- m/z sort order of output peaks
- 'keep_all' (non-random) sampling mode
- Empty MGF input
- Single-spectrum MGF
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'workflow'))
from generate_decoy import generate_decoy_mgf


# ---------------------------------------------------------------------------
# Minimal MGF fixtures
# ---------------------------------------------------------------------------

TWO_SPECTRUM_MGF = """\
BEGIN IONS
TITLE=spectrum1
PEPMASS=500.2345
CHARGE=2+
RTINSECONDS=120.5
100.00000 1000.00000
200.00000 2000.00000
300.00000 3000.00000
400.00000 4000.00000
END IONS
BEGIN IONS
TITLE=spectrum2
PEPMASS=612.4100
CHARGE=3+
RTINSECONDS=240.1
150.00000 1500.00000
250.00000 2500.00000
350.00000 3500.00000
450.00000 4500.00000
END IONS
"""

SINGLE_SPECTRUM_MGF = """\
BEGIN IONS
TITLE=single
PEPMASS=400.0
CHARGE=1+
100.00000 500.00000
200.00000 600.00000
END IONS
"""

EMPTY_MGF = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_mgf(path, content):
    with open(path, 'w') as f:
        f.write(content)


def parse_spectra(path):
    """Return list of dicts with 'meta' (list of str) and 'peaks' (list of (mz, intensity))."""
    spectra = []
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == 'BEGIN IONS':
                current = {'meta': [], 'peaks': []}
            elif line == 'END IONS':
                if current is not None:
                    spectra.append(current)
                current = None
            elif current is not None and line:
                if '=' in line:
                    current['meta'].append(line)
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        current['peaks'].append((float(parts[0]), float(parts[1])))
    return spectra


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestOutputFileCreation:
    def test_output_file_created(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        assert outfile.exists()

    def test_empty_mgf_creates_output(self, tmp_path):
        infile = tmp_path / "empty.mgf"
        outfile = tmp_path / "empty_decoy.mgf"
        write_mgf(infile, EMPTY_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        assert outfile.exists()
        assert parse_spectra(str(outfile)) == []


class TestSpectrumCount:
    def test_two_spectra_preserved(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        assert len(parse_spectra(str(outfile))) == 2

    def test_single_spectrum_preserved(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, SINGLE_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        assert len(parse_spectra(str(outfile))) == 1


class TestPeakCount:
    def test_peak_count_preserved_per_spectrum(self, tmp_path):
        """Decoy spectrum must contain the same total number of peaks as the original."""
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        for spectrum in spectra:
            assert len(spectrum['peaks']) == 4

    def test_single_spectrum_peak_count(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, SINGLE_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        assert len(spectra[0]['peaks']) == 2


class TestMetadata:
    def test_title_copied_to_output(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        assert any('TITLE=spectrum1' in m for m in spectra[0]['meta'])
        assert any('TITLE=spectrum2' in m for m in spectra[1]['meta'])

    def test_pepmass_copied(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        assert any('PEPMASS=500.2345' in m for m in spectra[0]['meta'])

    def test_charge_copied(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        assert any('CHARGE=2+' in m for m in spectra[0]['meta'])


class TestPeakOrder:
    def test_peaks_sorted_by_mz(self, tmp_path):
        """Output peaks must be sorted in ascending m/z order."""
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile))
        spectra = parse_spectra(str(outfile))
        for spectrum in spectra:
            mzs = [p[0] for p in spectrum['peaks']]
            assert mzs == sorted(mzs), f"Peaks not sorted: {mzs}"


class TestSamplingModes:
    def test_keep_all_mode_produces_output(self, tmp_path):
        """Non-random ('keep_all') mode: all original peaks retained, no noise injected."""
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile), peak_sampling='keep_all')
        spectra = parse_spectra(str(outfile))
        assert len(spectra) == 2

    def test_keep_all_mode_peak_count(self, tmp_path):
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile), peak_sampling='keep_all')
        spectra = parse_spectra(str(outfile))
        for spectrum in spectra:
            assert len(spectrum['peaks']) == 4

    def test_random_mode_default(self, tmp_path):
        """Default (random) mode should produce valid output."""
        infile = tmp_path / "sample.mgf"
        outfile = tmp_path / "sample_decoy.mgf"
        write_mgf(infile, TWO_SPECTRUM_MGF)
        generate_decoy_mgf(str(infile), str(outfile), sampling_rate=0.5, peak_sampling='random')
        spectra = parse_spectra(str(outfile))
        assert len(spectra) == 2
        for spectrum in spectra:
            assert len(spectrum['peaks']) == 4
