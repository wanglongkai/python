"""
    python类支持多继承
"""
class Student(object):
    __slots__ =  ('name', 'score', '__score', '__addr')
    age = 25 # 类属性，通过类名访问

    # 构造函数
    def __init__(self, name, score):
        self.name = name # 实例属性，通过实例访问
        self.score = score
        self.__score = score # 私有属性，只能类内部访问，外部不能访问

    '''
        getter setter
        只定义getter,就是一个只读属性
    '''
    @property
    def addr(self):
        return self.__addr

    @addr.setter
    def addr(self, value):
        self.__addr = value


"""
    枚举类
"""
from enum import Enum, unique

@unique # 保证去重，避免枚举重复
class Weekday(Enum):
    Sun = 0
    Mon = 1
    Tue = 2
    Wed = 3
    Thu = 4
    Fri = 5
    Sat = 6




if __name__ == '__main__':
    wlk = Student('wanglongkai', 100)

    print(isinstance(wlk, Student))
    wlk.addr = 'chongqing'
    print(wlk.addr)

    print(Weekday.Fri, Weekday.Fri.value)