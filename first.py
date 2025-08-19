class Cat(object):
  def say():
    print('i am a cat')

class Dog(object):
  def say():
    print('i am a dog')


animals = [Cat, Dog]
for animal in animals:
  animal.say()





