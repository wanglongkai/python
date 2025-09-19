# Python

## 环境搭建

1. 官网下载 python，点击安装，安装时记得勾选`add to path`;
2. 测试是否安装成功`py -V`, 打印版本号即成功

## 虚拟环境与激活

- 创建虚拟环境`py -3 -m venv .venv`
- 激活虚拟环境`.venv\Scripts\activate`

## pip 包写入 requirements.txt

`pip freeze > requirements.txt`

## 从 requirements.txt 安装

`pip install -r requirements.txt`

## uv 包管理工具

0. windows cmd 安装 uv

```shell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
set Path=C:\Users\admin\.local\bin;%Path%
```

1. 新建`pyproject.toml`

```toml
[project]
name = "projectName"
version = "0.0.1"
```

2. 添加新的依赖 `uv add fastapi`
3. 删除依赖`uv remove fastapi`
4. 如果是项目已有`pyproject.toml`, 可以通过`uv sync`命令安装所有依赖，同时会生成一个`uv.lock`文件
