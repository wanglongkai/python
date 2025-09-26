"""
    原始方案
"""
def originMethod():
    try:
        f = open('/path/file', 'r')
        print(f.read())
    finally:
        if f:
            f.close()


"""
    推荐方案： 不用手动处理文件报错及关闭文件
"""
def recommendMethod():
    with open('/path/file', 'r') as f:
        """
            f.read(): 一次读取所有内容
            f.read(size): 一次读取指定大小内容
            f.readline(): 每次读取一行
            f.readlines(): 读取所有内容，并以list形式返回
        """
        for line in f.readlines():
            print(line.strip()) # 把末尾的'\n'删掉


recommendMethod()