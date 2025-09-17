from fastapi import APIRouter, Form, File, UploadFile
from typing import Optional, List
from pydantic import BaseModel, Field
from datetime import date


userApi = APIRouter()

@userApi.post('/login/{name}', summary='用户登录接口 ')
def user_login(name, queryKey: str, queryDefault: Optional[str] = None):
    '''
        路径参数、query参数
    '''
    return {
        "name": name,
        "queryKey": queryKey,
        "queryDefault": queryDefault
    }

class Addr(BaseModel):
    province: str
    city: str

class User(BaseModel):
    name: str = Field(description='姓名')
    age: int = Field(default=0, lt=100, gt=0, description='年龄int类型，并且在0-100之间') 
    birth: Optional[date] = Field(default=None, description='生日')
    addr: Addr | None = Field(default=None, description='Addr 或者 None； 等同于 typing的Union')

@userApi.post('/bodydata', summary='请求体数据')
def test_body_data(user: User):
    '''
        body形式传送数据，并利用pydantic库进行数据格式校验
    '''
    return user

@userApi.post('/form/data', summary='表单提交')
def test_form_data(username: str = Form(), password: str = Form()):
    '''
        表单形式数据提交
    '''
    return {
        "username": username,
        "password": password
    }


@userApi.post('/smallfile', summary='适合小文件单文件上传')
def create_file(file: bytes = File()):
    return {
        'file_size': len(file)
    }

@userApi.post('/smallfiles', summary='适合小文件多文件上传')
def create_file(files: List[bytes] = File()):
    return {
        'fileslen': len(files)
    }


@userApi.post('/bigfile', summary='适合大文件单文件上传')
def create_file(file: UploadFile):
    return {
        'fileName': file.filename
    }

@userApi.post('/bigfiles', summary='适合大文件多文件上传')
def create_file(files: List[UploadFile]):
    return {
        'fileNames': [file.filename for file in files]
    }


