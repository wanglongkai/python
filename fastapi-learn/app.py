from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
import uvicorn
from userapi.userapi import userApi
from shopapi.shopapi import shopApi
from pdbapi.pdbapi import pdbApi

app = FastAPI(
    title="学习Python FastAPI项目",
    description="包含用户管理、商店管理和PDB分子转换功能的API服务",
    version="1.0.0"
)
# 设置静态资源文件夹
app.mount('/statics', StaticFiles(directory='statics'))

@app.get('/', tags=['根路由'], summary='根路由组件')
def read_root():
    return {'Hello': 'fastapi', 'message': '欢迎使用FastAPI学习项目'}

app.include_router(userApi, prefix='/user', tags=['用户管理接口'])
app.include_router(shopApi, prefix='/shop', tags=['商店管理接口'])
app.include_router(pdbApi, prefix='/pdb', tags=['PDB分子转换接口'])


if __name__ == '__main__':
    uvicorn.run('app:app', port=8080, reload=True)