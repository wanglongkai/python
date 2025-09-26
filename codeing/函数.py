"""
    位置参数
    默认参数
    关键字参数
    / 和 * 分隔三类参数 /之前的只能按位置传递， /和*之间的位置和关键字都可以， *之后的只能按关键字参数
"""

def testParams(name, age, /, addr, city='chongqing',color='red',*, hobit, cars=None):
    print(f'name={name}, age={age}, addr={addr}, city={city}, color={color}, hobit={hobit}, cars={cars}')




testParams('wlk',18,city='beijing',addr='haidian', color='blue', hobit=[1,2], cars=['car'])


def testKeyParams(name, age, *tuples, **keywords):
    print(name, age, tuples, keywords)


testKeyParams('wlk', 18, 'haha', 'dsfsdf', key1='keyw1', key2='keyw2')
testKeyParams('wlk',18, *[12,24], **{'key1': 'keys1', 'key2': 'keys2'}) # *展开列表、集合、元组，**展开字典