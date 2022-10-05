from this import d
import pandas as pd

def main():
    df = pd.read_csv("ds01_final.csv")
    df = df[[(("shelf"in str(de)) or ("clastic"in str(de)) ) for de in df["Depositional environment"]]]
    df = df[[((("meteoric" in em) and ("active" in em) and ("palaeo" in em))
            or (((("paeleo" in em)) and ("meteoric" in em)) or ("meteoric" in em)) and not ("active" in em))
            for em in df["Emplacement mechanism"]]]

    df = df[["Margin type", "Geology of aquifer", "Depositional environment", "Emplacement mechanism", "Reference"]]
    df.to_csv("subset.csv")
    df.describe()

if __name__=="__main__":
    main()