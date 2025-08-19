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



