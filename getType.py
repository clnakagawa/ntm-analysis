import numpy as np
import pandas as pd

def main():
    fname = "C:\\Users\cyn06\OneDrive - New York State Office of Information Technology Services\Downloads\lspn_gss_2023-06-20.csv"
    lpsn = pd.read_csv(fname)
    print(lpsn.head())

if __name__ == "__main__":
    main()