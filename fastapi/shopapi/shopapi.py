from fastapi import APIRouter

shopApi = APIRouter()

@shopApi.get('/shops', summary='获取购物清单信息')
def get_shops():
    return {'shop': 'shops'}


@shopApi.get('/shopInfo', summary='获取购物信息')
def get_shop_info():
    return {'info': 'info'}