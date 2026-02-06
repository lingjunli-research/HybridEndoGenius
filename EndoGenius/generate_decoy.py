import sys
import re
import numpy as np

def generate_decoy_mgf(input_mgf, output_mgf, sampling_rate=0.5, peak_sampling='random'):
    peaks_distr = []
    with open(input_mgf, 'r') as f:
        for line in f:
            line = line.strip()
            if line and line[0].isdigit():
                parts = re.split(r'\s+', line)
                if len(parts) >= 2:
                    try:
                        mz, intensity = map(float, parts[:2])
                        peaks_distr.append([mz, intensity])
                    except ValueError:
                        continue

    sampling_peaks_distr, noise_peaks_distr, decoy_peaks_distr = [], [], []

    with open(input_mgf, 'r') as f_in, open(output_mgf, 'w') as f_out:
        while True:
            line = f_in.readline()
            if not line:
                break
            if line.strip() == '':
                continue

            peak_list = []
            scan, peptide_mass = None, None
            while "END IONS" not in line:
                if 'BEGIN IONS' in line or '=' in line:
                    f_out.write(line)
                    line = f_in.readline()
                    continue

                line = line.strip()
                parts = re.split(r'\s+', line)
                if len(parts) >= 2:
                    try:
                        mz, intensity = map(float, parts[:2])
                        peak_list.append([mz, intensity])
                    except ValueError:
                        pass
                line = f_in.readline()

            num_peaks = len(peak_list)
            num_sampling = int(num_peaks * sampling_rate)
            num_noise = num_peaks - num_sampling
            np.random.seed(99)

            if peak_sampling == 'random':
                sampling_peaks = np.random.choice(len(peak_list), num_sampling, replace=False)
                sampling_peaks = [peak_list[i] for i in sampling_peaks]
                noise_peaks = np.random.choice(len(peaks_distr), num_noise, replace=False)
                noise_peaks = [peaks_distr[i] for i in noise_peaks]
            else:
                sampling_peaks = peak_list
                noise_peaks = []

            sampling_peaks_distr += sampling_peaks
            noise_peaks_distr += noise_peaks
            decoy_peaks_distr += sampling_peaks + noise_peaks

            sorted_peaks = sorted(sampling_peaks + noise_peaks, key=lambda x: x[0])
            for mz, intensity in sorted_peaks:
                f_out.write("{0:.5f} {1:.5f}\n".format(mz, intensity))
            f_out.write("END IONS\n\n")

if __name__ == "__main__":
    input_mgf = sys.argv[1]
    output_mgf = input_mgf.replace(".mgf", "_decoy.mgf")
    sampling_rate = 0.5
    peak_sampling = "random"
    generate_decoy_mgf(input_mgf, output_mgf, sampling_rate, peak_sampling)
