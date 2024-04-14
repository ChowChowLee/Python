import os
import re
import pandas as pd
import matplotlib.pyplot as plt

DATA_FOLDER = "C://Users//Lee//PycharmProjects//pythonProjectAlpacaBaby//Python//TPD//"

# Set your parameter of TPD plot
X_NAME = "'1/0'"
Y_NAME = "'0/0'"
X_LABEL = 'Temperature (K)'
Y_LABEL = 'Ion current (A)'
Y_MULTIPLE = 10.0

def show_file_list(DATA_FOLDER: str) -> list[str]:
    all_dta_file = [file for file in os.listdir(DATA_FOLDER) if file.endswith(".asc")]
    files_list = [i for i in all_dta_file]  # Show all the dta files there
    print(files_list)

FILE_NAME_LIST = [
                'coverage dependence of o2  ramp4 exp1.asc',
                'coverage dependence of o2  ramp4 exp2.asc',
                'coverage dependence of o2  ramp4 exp3.asc',
                'coverage dependence of o2  ramp4 exp4.asc'
                ]

LEGEND_LIST = [
                 'exp1',
                 'exp2',
                 'exp2',
                 'exp3',
              ]

HEADER_LINE_HEAD = 'Cycle' + '\t'

def load_data(file_path: str, header_line_head: str) -> pd.DataFrame:
    with open(file_path, 'r') as file:
        while True:
            line = file.readline()
            if line[:len(header_line_head)] == header_line_head:
                header = line.strip().split('\t')
                break
        text = [line.strip().split('\t') for line in file.readlines() if re.search(r'[0-9]', line[0])]
    return pd.DataFrame(text, columns=header)


def plot_spectrum(FILE_NAME_LIST: list[str], LEGEND_LIST: list[str], DATA_FOLDER: str, HEADER_LINE_HEAD: str,
                  X_NAME: str, Y_NAME: str, Y_MULTIPLE: float) -> None:
    plt.figure(figsize=(8, 6))

    for file_name, i_legend in zip(FILE_NAME_LIST, LEGEND_LIST):
        try:
            dt = load_data(DATA_FOLDER + file_name, HEADER_LINE_HEAD)
            plt.plot(dt[X_NAME].apply(lambda x: float(x)), dt[Y_NAME].apply(lambda x: float(x)) * Y_MULTIPLE,
                     label=i_legend)
            print(dt.head())
        except:
            print('No column name: ' + X_NAME + ' or ' + Y_NAME + ' in ' + file_name + '\nonly :' + ', '.join(
                dt.columns))

    plt.xlabel(X_LABEL)
    plt.ylabel(Y_LABEL)
    plt.xlim(left=300, right=1380)
    plt.ylim(bottom=1e-9, top=4.0e-9)
    plt.legend(loc='upper left', bbox_to_anchor=(0.8, 0.9))
    plt.title("TPD")
    plt.show()


def main() -> None:
    show_file_list(DATA_FOLDER)
    plot_spectrum(FILE_NAME_LIST, LEGEND_LIST, DATA_FOLDER, HEADER_LINE_HEAD, X_NAME, Y_NAME, Y_MULTIPLE)

if __name__ == "__main__":
    main()