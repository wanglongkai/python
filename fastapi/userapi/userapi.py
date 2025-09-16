from fastapi import APIRouter

userApi = APIRouter()

@userApi.post('/login', summary='用户登录接口')
def user_login():
    return {'user': 'login'}