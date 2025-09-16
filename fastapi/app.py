from fastapi import FastAPI
import uvicorn
from userapi.userapi import userApi
from shopapi.shopapi import shopApi

app = FastAPI()

@app.get('/', tags=['根路由'], summary='根路由组件')
def read_root():
    return {'Hello': 'fastapi'}

app.include_router(userApi, prefix='/user', tags=['这是user的接口'])
app.include_router(shopApi, prefix='/shop', tags=['这是shop的接口'])


if __name__ == '__main__':
    uvicorn.run('app:app', port=8080, reload=True)