#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 清理Raser更新前后文件夹
rule：
    * 清理不以'_out'的所有文件夹（Raser新版结果文件夹都以'_out'结尾）
    * 清理以'_out'结尾的空文件夹
"""

import os
import sys
import shutil

def main():
    print("Usage: python3 clean_folder.py [path] \n{example}\n".format(
        example="Example: python3 clean_folder.py ./PRJEB20313")
    )
    abs_path = os.path.abspath(os.path.expanduser(sys.argv[1]))
    for path in os.listdir(abs_path):
        path = os.path.join(abs_path, path)
        if os.path.isfile(path):
            continue
        for folder in os.listdir(path):
            if os.path.isfile(os.path.join(path, folder)):
                continue
            if not folder.endswith("_out"):
                print("clean dir: {0}".format(os.path.join(path, folder)))
                shutil.rmtree(os.path.join(path, folder))




if __name__ == '__main__':
    main()
