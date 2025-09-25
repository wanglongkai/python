# 数据类型知识点

- python 默认是 Unicode 编码
  str 有编解码方法：`encode('utf-8'), decode('utf-8')`  
  比如：  
  `b'\xe4\xb8\xad\xe6\x96\x87'.decode('utf-8')`结果就是'中文'。（bytes）

- f-string  
  f('你好{name}, 我身高{hight:.2f}')
