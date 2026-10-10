import re
import settings

periodic_table = settings.periodic_table


def parse_formula(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*\.\d+|\d*)', formula)
    element_dict = {}
    for element, amount in elements:
        amount = float(amount) if amount else 1.0
        element_dict[element] = amount
    return element_dict


def wt_to_at_batch(formulas):
    results_list = []

    for input_formula in formulas:
        input_formula_dict = parse_formula(input_formula)

        result_dict = {}
        result_dict['formula'] = input_formula
        result_dict['elements'] = input_formula_dict

        n_list = [input_formula_dict[element] / periodic_table.get(element, 1) for element in input_formula_dict]
        n_all = sum(n_list)

        atom_and_at_percent_dict = {element: n / n_all * 100 for element, n in zip(input_formula_dict, n_list)}
        result_dict['at_percent'] = atom_and_at_percent_dict
        result_dict['input_sum'] = sum(input_formula_dict.values())
        result_dict['output_sum'] = sum(atom_and_at_percent_dict.values())

        results_list.append(result_dict)

    return results_list


if __name__ == "__main__":
    # 在这里输入你的化学式列表
    formulas = [
        "Ga61In25Sn13Zn1",
        "Ga68.5In21.5Sn10",
        "In51Bi32.5Sn16.5",
        "Bi58Sn42",
        "In25.2Sn17.3Bi57.5",
        "In51.34Sn15.56Bi33.1",
        "In4Sn40Bi56",
        "Sn22Bi50Pb28",
        "Sn15.5Bi52.5Pb32",
        "Sn26Bi53Cd21",
        "Bi51.6Cd8.2Pb40.2",
        "Sn51.2Cd30.6Pb18.2",
        "In21Sn12Bi49Pb18",
        "Sn13.3Bi50Cd10Pb26.7",
        "Bi44.7Pb22.6In19.1Sn8.3Cd5.3",
        "Bi49In21Pb18Sn12",
        "Bi31.6In48.8Sn19.6",
        "Bi50Pb26.7Sn13.3Cd10",
        "Bi52Pb30Sn18",
        "Sn91Zn9",
        "Ga67In20.5Sn12.5",
        "Ga86.3Sn10.8Zn2.9",
        "Ga86.08425Sn10.773Zn2.89275Pb0.25",
        "Ga85.8685Sn10.746Zn2.8855Pb0.5",
        "Ga85.437Sn10.692Zn2.871Pb1",
        "Ga84.574Sn10.584Zn2.842Pb2",
        "Ga81.985Sn10.26Zn2.755Pb5",
        "Ga86.29137Sn10.79892Zn2.89971Li0.01",
        "Ga86.25685Sn10.7946Zn2.89855Li0.05",
        "Ga86.2137Sn10.7892Zn2.8971Li0.1",
        "Ga85.8685Sn10.746Zn2.8855Al0.5",
        "Ga85.437Sn10.692Zn2.871Al1",
        "Ga84.574Sn10.584Zn2.842Al2",
        "Ga85.437Sn10.692Zn2.871In1",
        "Ga81.985Sn10.26Zn2.755In5",
        "Ga77.67Sn9.72Zn2.61In10",
        "Ga86.29137Sn10.79892Zn2.89971Na0.01",
        "Ga86.25685Sn10.7946Zn2.89855Na0.05",
        "Ga86.2137Sn10.7892Zn2.8971Na0.1",
        "Ga85.8685Sn10.746Zn2.8855Na0.5",
        "Ga85.437Sn10.692Zn2.871Na1",
        "Ga60.7In25.2Sn13Zn1.1",
        "Ga78In22",
        "Ga86.5Sn13.5",
        "Bi60Cd40",
        "Bi57Sn43",
        "Bi49In51",
        "Bi55.5Pb44.5",
        "Cd25.6In74.4",
        "Ga75.5In24.5",
        "Ga67.8318In21.4Sn8.4888Zn2.2794",
        "Ga78.3In14.9Sn6.8",
        "Ga68.363In21.457Sn9.98Gd0.2",
        "Ga67.678In21.242Sn9.88Gd1.2",
        "Ga66.9245In21.0055Sn9.77Gd2.3",
        "Ga90Zn10",
        "Ga62Sn32Bi6",
        "Ga91.6Sn8.4",
        "Ga90In10",
        "Ga80In20",
        "Ga70In30",
        "Ga60In40",
        "Ga50In50",
        "Ga40In60",
        "Ga30In70",
        "Ga90Sn10",
        "Ga80Sn20",
        "Ga70Sn30",
        "Ga60Sn40",
        "Ga50Sn50",
        "Ga80Zn20",
        "Ga70Zn30",
        "Ga80In5Sn15",
        "Ga80In10Sn10",
        "Ga80In15Sn5",
        "Ga80In7.5Sn12.5",
        "Ga80In12.5Sn7.5",
        "Ga77.6In9.7Sn9.7Zn3",
        "Ga77.6In9.7Sn9.7Ge3",
        "Ga77.6In9.7Sn9.7Al3",
        "Ga68In20Sn12",
        "Bi45Pb23In19Sn8Cd5",
        "Bi31.6In48.8Sn19.6",
        "Cd32.3Sn67.7",
        "Pb37Sn63",
        "Sn57Tl43",
        "Ga86Sn11Zn3",
        "Sn38Bi50In12",
        "Sn41Bi47In12",
        "Sn47Bi41In12",
        "Sn44Bi44In12",
        "Sn50Bi38In12",
        "Sn53Bi35In12",
        "Ga86.3Sn10.8Zn2.9",
        "Ga86.08425Sn10.773Zn2.89275Bi0.25",
        "Ga85.8685Sn10.746Zn2.8855Bi0.5",
        "Ga85.437Sn10.692Zn2.871Bi1",
        "Ga84.574Sn10.584Zn2.842Bi2",
        "Ga81.985Sn10.26Zn2.755Bi5",
        "Sn41.95Bi58Ti0.05",
        "Sn41.9Bi58Ti0.1",
        "Sn41.8Bi58Ti0.2",
        "Sn41.7Bi58Ti0.3",
        "Sn41.5Bi58Ti0.5",
        "Sn41Bi58Ti1",
        "Sn40Bi58In2",
        "Sn40Bi58Ag2",
        "Sn38Bi58In2Ag2",
        "Sn41.5Bi58La0.5",
        "Sn41.5Bi58Ag0.5",
        "Sn58Bi40Zn2",
        "Sn57.995Bi40Zn2Al0.005",
        "Sn59.9Bi40Cu0.1",
        "Sn57.9Bi40Cu0.1Zn2",
        "Sn64.7Bi35Ag0.3",
        "Sn65Bi35Ag1",
        "Sn41Bi58In1",
        "Sn40Bi58In2",
        "Sn39Bi58In3",
        "Sn38Bi58In4",
        "Sn37Bi58In5",
        "Sn46.2Bi52Sb1.8",
        "Sn50.2Bi48Sb1.8",
        "Sn54.2Bi44Sb1.8",
        "Sn51Bi48Sb1.0",
        "Sn50Bi48Sb2",
        "Sn49.2Bi48Sb2.8",
        "Sn41.9Bi58Ag0.1",
        "Sn42Bi57.6Ag0.4",
        "Sn41.5Bi58Ag0.5",
        "Sn41Bi58Ag1",
        "Sn40Bi58Ag1In1",
        "Sn38Bi58Ag1In3",
        "Sn38Bi58Ag4",
        "Sn35Bi58In3Ag4",
        "Sn42Bi58Co0.05",
        "Sn42Bi58Co0.5",
        "Sn59.9Bi40Cu0.1",
        "Sn57.9Bi40Zn2Cu0.1",
        "Sn40.7Bi56.2Sb3.1",
        "Sn41Bi53.2Sb5.8",
        "Sn50.6Bi48Sb1.4",
        "Sn49.6Bi48Sb2.4",
        "Sn70Bi20In10",
    ]

    results = wt_to_at_batch(formulas)
    list1 = []
    # 打印结果列表
    for result in results:
        # print(result['at_percent'])
        list1.append(result['at_percent'])
    print(list1)

    import pandas as pd

    # Your list of dictionaries
    data = list1

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Add the chemical formula column
    df['Chemical Formula'] = [
        ''.join([f'{element}{amount}' for element, amount in row.items() if element != 'Chemical Formula']) for row in
        data]

    df = df.fillna(0)

    # Reorder columns
    columns_order = ['Chemical Formula'] + list(df.columns[:-1])
    df = df[columns_order]

    # Save to Excel
    df.to_excel(r'C:\Users\senji\Desktop\AT%.xlsx', index=False)


