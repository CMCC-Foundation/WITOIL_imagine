import argparse
import pandas as pd
import os


def read_oilbase(option, value, data_path, output_path) -> str:

    data = pd.read_csv(data_path, sep=",", comment="#")

    if option == "name":
        try:  # collecting line with correspondent oil name
            data.index = data.name
            oil = data.loc[value, :]

        except:  # if there is no exact match
            # line below gets all names that start with the same letter as input
            oil_withletter = data[
                data["name"].str.startswith(value[0])
            ].name.values.tolist()

            print("There is no oil with that name. Do you mean:")
            print("\n".join(oil_withletter))

            raise ValueError(
                "No oil with that name in database. Try suggestions mentioned above."
            )

    else:  # uses closest API value to return an oil
        data.index = data.name
        oil = data.loc[(data.api - float(value)).abs().idxmin()]

    # writing oil values to oil_file.txt that will be use in medslik-II simulation
    output_file = os.path.join(output_path, "oil_file.txt")
    with open(output_file, "w") as f:
        f.write(f"{oil.name}\n")
        f.write(f"{oil.api:.02f}          API of Oil\n")
        f.write(f"{oil.density:.03f}          Density of Oil\n")
        f.write(f"{oil.residual_density:.03f}          Residual Density of Oil\n")
        f.write(f"{oil.residual_percentage:.02f}          Residual Percent of Oil\n")
        f.write(f"{oil.viscosity:.02f}          Viscosity of Oil\n")
        f.write(
            f"{oil.temperature:.02f}          Temperature at which Viscosity determined\n"
        )
        f.write(f"{oil.vapour_pressure:.03f}          Vapour Pressure of Oil (bar)\n")

        f.close()
    return output_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Collect oil information from the oilbase"
    )
    parser.add_argument("option", type=str, help='Either "name" or "api"')
    parser.add_argument(
        "value", help="Oil name or Oil API. Oil name must be exact value"
    )
    parser.add_argument("output_folder", help="path to experiment output folder")
    args = parser.parse_args()
    home_path = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    filename = "oilbase.csv"
    oil_path = os.path.join(home_path, "data", filename)
    print(f"Reading oil from database {oil_path}")
    read_oilbase(args.option, args.value, oil_path, args.output_folder)
