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
        'Sn62Bi32Ga6',
    ]

    results = wt_to_at_batch(formulas)

    # 打印结果列表
    for result in results:
        print(f'化学式：{result["formula"]}')
        print(f'你输入的元素成分为：{result["elements"]}')
        print(f'转化成at%为：{result["at_percent"]}')
        print(f'输入字典值的求和：{result["input_sum"]}')
        print(f'输出字典值的求和：{result["output_sum"]}')
        print('\n')
