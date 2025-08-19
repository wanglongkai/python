# Python
## 环境搭建
1. 官网下载python，点击安装，安装时记得勾选`add to path`;
2. 测试是否安装成功`py -V`, 打印版本号即成功

## 虚拟环境与激活
- 创建虚拟环境`py -3 -m venv .venv`
- 激活虚拟环境`.venv\Scripts\activate`

## pip包写入requirements.txt
`pip freeze > requirements.txt`

## 从 requirements.txt 安装
`pip install -r requirements.txt`

## 常见数据类型
- None 全局唯一`id()函数比较两个None是相等的`
- int、float、str、list


## 流程控制
- if、elif、else
- for
- while
- break、continue
- pass
- match 类似与javascript中的switch
