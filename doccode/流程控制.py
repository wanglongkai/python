"""
流程控制
"""

# 条件语句
x = 1
if x == 1:
    print('x is 1')
elif x == 2:
    print('x is 2')
else:
    print('x is not 1 or 2')


# 循环语句
users = {
	'admin': '123456',
	'user1': '123456',
}

for username, password in users.items():
    print(username, password)

for username, password in users.copy().items():
    print(username, password)
    if username == 'admin':
        users[username] = '123456789'
print(users)

for x in range(10):
    print(x)
    if x == 5:
        break
else:
    print('for loop is done')


print(type({1,2,3}))


# 函数参数
def f(pos1, pos2, /, pos_or_kwd, *, kwd1, kwd2):
    print(pos1, pos2, pos_or_kwd, kwd1, kwd2)

# 调用函数
f(1, 2, 'a', kwd1=100, kwd2=200)






