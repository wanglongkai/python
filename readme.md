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

## 常见数据类型

- 列表 [1,2,3]
- 元组 (1,2,3)
- 集合 {1,2,3}
- 字典 {'a':1, 'b':2}
- 字符串 'hello'
- 布尔值 True False
- 空值 None
- 数字 1 2 3 4 5 6 7 8 9 0
- 浮点数 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 0.0

## 抽象基类

```python
# 定义抽象基类
class Base:
  def say(self):
    # 该方法不实现，会报错
    raise NotImplementedError

  def hello(self):
    # 该方法可以不实现
    pass

# 实现抽象基类
class Cat(Base):
  pass

cat = Cat()
cat.say()
cat.hello()
```

`isinstance()、 type()`

## 循环的 else 子句

如果循环在未执行 break 的情况下结束，会执行 else 子句；
如果循环执行了 break，会跳过 else 子句。

```python
for x in range(10):
    print(x)
    if x == 5:
        break
else:
    print('for loop is done')
```

## 解包

\* 解包列表或元祖
\*\* 解包字典
