import os
import re
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional, Union

DATA_FOLDER = ("C://Users//Lee//PycharmProjects//pythonProjectAlpacaBaby//Python//XRD//")
HEADER_LINE_HEAD = None

# Set your parameter
X_NAME = 0
Y_NAME = 1
X_LABEL = '2θ'
Y_LABEL = 'Relative intensity'
YMAX_MULTIPLIER = 5
Y_MULTIPLE =1.0

def show_file_list(DATA_FOLDER: str) -> None:
    all_dta_file = [file for file in os.listdir(DATA_FOLDER) if file.endswith(".xy")]
    [print(i+'\n') for i in all_dta_file]  # Show all the dta files there

FILE_NAME_LIST = [
                'OER Sample 1 (Ti TiO2 Sb2S3)_25112020_bg_subtracted.xy',
                'OER Sample 2 (Ti TiO2 Sb2S3 MoO3)_25112020_bg_subtracted.xy',
                'OER Sample 3(2) (Ti TiO2 Sb2S3 Ni)_25112020_bg_subtracted.xy',
                'OER Sample 4 (Ti TiO2 Sb2S3 MoO3 Ni)_25112020_bg_subtracted.xy'
                ]

LEGEND_LIST = [
                'OER Sample 1 (Ti/TiO2/Sb2S3)',
                'OER Sample 2 (Ti/TiO2/Sb2S3/MoO3)',
                'OER Sample 3 (Ti/TiO2/Sb2S3/Ni)',
                'OER Sample 4 (Ti/TiO2/Sb2S3/MoO3/Ni)'
              ]

def lr_strip_split(line: str) -> str:
    return line.rstrip('\n').split(' ')

def load_data(file_path: str, header_line_head: Optional[str] = HEADER_LINE_HEAD) -> pd.DataFrame:
    with open(file_path, 'r') as file:
        if header_line_head is None:
            text = core_logic(file, header_line_head)
            dt = pd.DataFrame(text)
            return dt
        text, header = core_logic(file, header_line_head)
        dt = pd.DataFrame(text, columns=header)
    return dt

def core_logic(file, header_line_head: Optional[str] = HEADER_LINE_HEAD) -> str:
    line = file.readline()
    if header_line_head is None:
        text = [lr_strip_split(line) for line in file.readlines() if re.search(r'[0-9]', line[:2])]
        return text
    if line[:len(header_line_head)] == header_line_head:
        print(line)
        header = lr_strip_split(line)
    text = [lr_strip_split(line) for line in file.readlines() if re.search(r'[0-9]', line[:2])]
    return text, header

def plot_spectrum(FILE_NAME_LIST: list[str],
                  LEGEND_LIST: list[str],
                  DATA_FOLDER: str,
                  HEADER_LINE_HEAD: Optional[str],
                  X_NAME: Optional[Union[str, int]],
                  Y_NAME: Optional[Union[str, int]],
                  Y_MULTIPLE: float) -> None:

    plt.figure(figsize=(6,4))

    for file_name, i_legend in zip(FILE_NAME_LIST, LEGEND_LIST):
        try:
            dt = load_data(DATA_FOLDER + file_name, HEADER_LINE_HEAD)
            ymax_previous = YMAX_MULTIPLIER * max(dt[1].astype(float))
            plt.plot(dt[X_NAME].astype(float), dt[Y_NAME].astype(float) * Y_MULTIPLE + ymax_previous, label=i_legend)
            ymax_previous += YMAX_MULTIPLIER * max(dt[1].astype(float))
        except:
            print('No column name: ' + X_NAME + ' or ' + Y_NAME + ' in ' + file_name + '\nonly :' + ', '.join(dt.columns))

    plt.xlabel(X_LABEL)
    plt.ylabel(Y_LABEL)
    plt.tick_params(labelleft=False, left=False)
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, 0.9))
    plt.title("XRD")
    plt.show()


def main() -> None:
    show_file_list(DATA_FOLDER)
    plot_spectrum(FILE_NAME_LIST, LEGEND_LIST, DATA_FOLDER, HEADER_LINE_HEAD, X_NAME,Y_NAME,Y_MULTIPLE)

if __name__ == "__main__":
    main()