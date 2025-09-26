def foo(s):
    return 10 / int(s)

def bar(s):
    return foo(s) * 2

def main():
    try:
        bar('0')
    except Exception as e:
        print('Error', e)
    finally:
        print('finally...')

    print('继续执行后续代码')


main()