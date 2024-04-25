import os
import re
import pandas as pd
import matplotlib.pyplot as plt

DATA_FOLDER = "C://Users//Lee//PycharmProjects//pythonProjectAlpacaBaby//Python//CyclicVoltammetry//"
HEADER_LINE_HEAD = '\t' + 'Pt'

# Set your parameter
X_NAME = 'Vf'
Y_NAME = 'Im'
X_LABEL = 'E / V vs. Ag / Agcl'  # 'Time (s)'
Y_LABEL = 'I / µA'  # 'Potential (mV)'
Y_MULTIPLE = 1000000

def show_file_list(DATA_FOLDER: str) -> None:
    all_dta_file = [file for file in os.listdir(DATA_FOLDER) if file.endswith(".DTA") or file.endswith("dta")]
    [print(i+'\n') for i in all_dta_file]  # Show all the dta files there

FILE_NAME_LIST = [
                'Sb2S3 MoO3 Ni buffer+urea TWLCV 0-1.5V 50mVs .DTA',
                'Sb2S3 Ni buffer+urea light TWLCV 0-1.5V 50mVs .DTA',
                'Sb2S3 Ni buffer+urea TWLCV 0-1.5V 50mVs .DTA'
                ]

LEGEND_LIST = [
                'Sb2S3 MoO3 Ni buffer+urea',
                'Sb2S3 Ni buffer+urea light',
                'Sb2S3 Ni buffer+urea'
              ]

def lr_strip_split_line(line: str) -> str:
    return line.rstrip('\n').lstrip('\t').split('\t')

def load_data(file_path: str, header_line_head: str = HEADER_LINE_HEAD) -> pd.DataFrame:
    with open(file_path, 'r') as file:
        text, header = core_logic(file, header_line_head)
    dt = pd.DataFrame(text, columns=header)
    return dt

def core_logic(file, header_line_head: str = HEADER_LINE_HEAD) -> str:
    while True:
        line = file.readline()
        if line[:len(header_line_head)] == header_line_head:
            print(line)
            header = lr_strip_split_line(line)
            break
    text = [lr_strip_split_line(line) for line in file.readlines() if re.search(r'\t[0-9]', line[:2])]
    return text, header

def fix_comma_to_float(x: str) -> float:
    return float(x.replace(',','.'))

def plot_spectrum(FILE_NAME_LIST: list[str],
                  LEGEND_LIST: list[str],
                  DATA_FOLDER: str,
                  HEADER_LINE_HEAD: str,
                  X_NAME: str,
                  Y_NAME: str,
                  Y_MULTIPLE: float) -> None:

    plt.figure(figsize=(6,4))

    for file_name, i_legend in zip(FILE_NAME_LIST, LEGEND_LIST):
        try:
            dt = load_data(DATA_FOLDER + file_name, HEADER_LINE_HEAD)
            plt.plot(dt[X_NAME].apply(fix_comma_to_float), dt[Y_NAME].apply(fix_comma_to_float) * Y_MULTIPLE, label=i_legend)
        except:
            print('No column name: ' + X_NAME + ' or ' + Y_NAME + ' in ' + file_name + '\nonly :' + ', '.join(dt.columns))

    plt.xlabel(X_LABEL)
    plt.ylabel(Y_LABEL)
    plt.legend(loc='upper left', bbox_to_anchor=(0.0,0.9))
    plt.title("Cyclic Voltammetry")
    plt.show()


def main() -> None:
    show_file_list(DATA_FOLDER)
    plot_spectrum(FILE_NAME_LIST, LEGEND_LIST, DATA_FOLDER, HEADER_LINE_HEAD, X_NAME, Y_NAME, Y_MULTIPLE)

if __name__ == "__main__":
    main()