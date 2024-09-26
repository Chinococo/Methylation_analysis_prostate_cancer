import csv
import os
from itertools import combinations

import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt



class Predict:
    def __init__(self,test_name:str):
        self.targetDir = "Predict/"+test_name
        self.choose = 4
        self.draw = []
        if os.path.exists(f"{self.targetDir}/result.csv"):
            os.remove(f"{self.targetDir}/result.csv")
        self.target = []
        self.testfile = pd.read_csv(f"{self.targetDir}/all_beta_normalized.csv", header='infer', index_col='Unnamed: 0')
        self.ansfile = pd.read_csv(f"{self.targetDir}/ans.csv")
        self.AnsTumor = set(self.ansfile[self.ansfile["state"] == "Tumor"]["id"].astype(str))
        self.AnsNormal = set(self.ansfile[self.ansfile["state"] == "Normal"]["id"].astype(str))

    def plot_data(self, data, threshold, labels, name, targetDir):
        fig, ax = plt.subplots()

        # 设置变量来跟踪是否已经添加图例条目
        added_labels = {'red': False, 'green': False}

        # 遍历第一行的 Series 数据
        for idx, (label, value) in enumerate(data.iloc[0].items()):  # 使用 items() 而不是 iteritems()
            if label in self.AnsTumor:
                # 仅在首次添加时为散点图提供图例
                if not added_labels['red']:
                    ax.scatter(idx, value, color='red', label='Tumor')
                    added_labels['red'] = True
                else:
                    ax.scatter(idx, value, color='red')
            else:
                if not added_labels['green']:
                    ax.scatter(idx, value, color='green', label='Non-Tumor')
                    added_labels['green'] = True
                else:
                    ax.scatter(idx, value, color='green')

        # 添加阈值线
        ax.axhline(y=threshold, color='black', linestyle='--', label='Threshold')

        # 设置 x 和 y 标签
        ax.set_xlabel('Index')
        ax.set_ylabel('Data')

        # 设置图标题
        ax.set_title(f'{name} in {targetDir.split("_")[-1]}')

        # 显示图例，避免重复条目
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # 显示网格
        plt.grid(True)

        # 展示图表
        plt.show()
    def LoadSelect(self):
        self.target.append("Predict/select_target.csv")


    def WriteToCsv(self, best_combination, best_acc, best_error, f1):
        with open('results.csv', 'a') as file:
            file.write(f"{best_combination},{best_acc:.2f}%,{best_error},{f1:.2f}\n")

    def Predict(self):
        for testpath in self.target:
            # Initialize temporary variables
            judgefile = pd.read_csv(testpath)
            cpg_list = judgefile["CpG"].tolist()

            # Generate all possible combinations of CpG sites
            all_combinations = []
            for i in range(1, self.choose + 1):
                all_combinations.extend(combinations(cpg_list, i))

            best_acc = 0
            best_combination = []
            best_error = []
            best_f1 = 0

            for selctCpg in all_combinations:
                if len(selctCpg) != self.choose:
                    continue
                print(selctCpg)
                Tumor = []
                for cpg in selctCpg:
                    TargetCutPoint = float(judgefile[judgefile["CpG"] == cpg]["cutpoint"].iloc[0])
                    test = self.testfile.loc[self.testfile.index == cpg]
                    if len(test) == 0:
                        print("不存在", cpg)
                        continue
                    selected_columns = test.columns[(test > TargetCutPoint).any()]
                    Tumor.extend(list(selected_columns))
                    if cpg not in self.draw:
                        self.draw.append(cpg)
                        self.plot_data(test, TargetCutPoint, test.columns, cpg, self.targetDir)

                element_counts = Counter(Tumor)
                result = [element for element, count in element_counts.items() if count > (self.choose / 2)]
                acc = len(self.AnsTumor.intersection(set(result))) / len(self.AnsTumor) * 100
                error = sorted(list(self.AnsTumor - self.AnsTumor.intersection(set(result))))

                # Calculate precision, recall, and F1 score
                tp = len(self.AnsTumor.intersection(set(result)))
                fp = len(set(result) - self.AnsTumor)
                fn = len(self.AnsTumor - set(result))

                if tp + fp == 0:
                    precision = 0
                else:
                    precision = tp / (tp + fp)

                if tp + fn == 0:
                    recall = 0
                else:
                    recall = tp / (tp + fn)

                if precision + recall == 0:
                    f1 = 0
                else:
                    f1 = 2 * (precision * recall) / (precision + recall)

                if acc > best_acc:
                    best_acc = acc
                    best_combination = selctCpg
                    best_error = error
                    best_f1 = f1

            self.WriteToCsv(best_combination, best_acc, best_error, best_f1)
            error_acc = (len(self.AnsNormal) - len(list(self.AnsNormal - set(result)))) / len(self.AnsNormal) * 100
            print(f"{best_combination},{best_acc:.2f}%,{best_error},{error_acc:.2f}%,{best_f1:.2f}")
            print(self.AnsNormal - set(self.AnsNormal - set(result)))
