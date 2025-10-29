from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import JSONResponse
from typing import Optional
from pydantic import BaseModel, Field
import tempfile
import os
import re
import rdkit
from rdkit import Chem

pdbApi = APIRouter()

class PDBResponse(BaseModel):
    """PDB转SMILES响应模型"""
    success: bool = Field(description="处理是否成功")
    smiles: Optional[str] = Field(default=None, description="生成的SMILES字符串")
    pocket_info: Optional[str] = Field(default=None, description="使用的口袋信息")
    message: str = Field(description="处理结果消息")
    temp_file: Optional[str] = Field(default=None, description="生成的临时文件路径")

def parse_pocket_info(pocket_str: str):
    """
    解析口袋信息字符串
    格式: A:1E8:701 (链:配体名:残基号)
    返回: (chain, ligand_name, residue_number)
    """
    pattern = r'^([A-Z]):([A-Z0-9]+):(\d+)$'
    match = re.match(pattern, pocket_str)
    if match:
        return match.group(1), match.group(2), int(match.group(3))
    else:
        raise ValueError(f"口袋信息格式错误: {pocket_str}，正确格式应为 链:配体名:残基号 (如: A:1E8:701)")

def extract_pocket_from_pdb_content(pdb_content: str, pocket_info: str):
    """
    从PDB内容中提取指定口袋区域的内容
    
    参数:
    pdb_content: PDB文件内容字符串
    pocket_info: 口袋信息字符串，格式为 "A:1E8:701"
    
    返回:
    过滤后的PDB内容字符串
    """
    # 解析口袋信息
    chain, ligand_name, residue_number = parse_pocket_info(pocket_info)
    
    # 过滤出指定口袋区域的HETATM记录
    filtered_lines = []
    found_hetatm = False
    
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            # 解析HETATM行的格式
            if len(line) >= 26:
                line_chain = line[21:22].strip()
                line_residue = line[17:20].strip()
                line_residue_num = line[22:26].strip()
                
                # 检查是否匹配指定的口袋信息
                if (line_chain == chain and 
                    line_residue == ligand_name and 
                    line_residue_num == str(residue_number)):
                    filtered_lines.append(line)
                    found_hetatm = True
        elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK', 'CRYST1', 'ORIGX', 'SCALE', 'MASTER', 'END')):
            # 保留PDB文件头信息
            filtered_lines.append(line)
    
    if not found_hetatm:
        raise ValueError(f"未找到匹配的口袋信息: {pocket_info}")
    
    return '\n'.join(filtered_lines)

def process_pdb_to_smiles(pdb_content: str, pocket_info: Optional[str] = None):
    """
    处理PDB内容并转换为SMILES
    
    参数:
    pdb_content: PDB文件内容
    pocket_info: 可选的口袋信息
    
    返回:
    (smiles_string, temp_file_path, message)
    """
    try:
        # 创建临时文件
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as temp_file:
            if pocket_info:
                # 使用口袋信息过滤
                try:
                    filtered_content = extract_pocket_from_pdb_content(pdb_content, pocket_info)
                    temp_file.write(filtered_content)
                    temp_file_path = temp_file.name
                    message = f"成功提取口袋信息 {pocket_info} 的配体分子"
                except ValueError as e:
                    # 如果口袋信息不匹配，使用完整文件
                    temp_file.write(pdb_content)
                    temp_file_path = temp_file.name
                    message = f"警告: {str(e)}，使用完整PDB文件进行处理"
            else:
                # 使用完整PDB文件
                temp_file.write(pdb_content)
                temp_file_path = temp_file.name
                message = "使用完整PDB文件进行处理"
        
        # 使用RDKit处理PDB文件
        mol = Chem.rdmolfiles.MolFromPDBFile(temp_file_path)
        
        if mol is None:
            raise ValueError("无法从PDB文件创建分子对象，请检查文件格式")
        
        # 转换为SMILES
        smiles = Chem.MolToSmiles(mol)
        
        return smiles, temp_file_path, message
        
    except Exception as e:
        # 清理临时文件
        if 'temp_file_path' in locals() and os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
        raise e

@pdbApi.post('/convert', response_model=PDBResponse, summary='PDB文件转SMILES')
async def convert_pdb_to_smiles(
    file: UploadFile = File(..., description="PDB文件"),
    pocket_info: Optional[str] = Form(None, description="口袋信息，格式: A:1E8:701 (链:配体名:残基号)")
):
    """
    将PDB文件转换为SMILES字符串
    
    参数:
    - **file**: PDB文件 (必需)
    - **pocket_info**: 口袋信息，格式为 "A:1E8:701" (可选)
      - A: 链标识符
      - 1E8: 配体名称
      - 701: 残基号
    
    返回:
    - **success**: 处理是否成功
    - **smiles**: 生成的SMILES字符串
    - **pocket_info**: 使用的口袋信息
    - **message**: 处理结果消息
    - **temp_file**: 生成的临时文件路径（用于调试）
    """
    temp_file_path = None
    
    try:
        # 验证文件类型
        if not file.filename.lower().endswith('.pdb'):
            raise HTTPException(status_code=400, detail="文件必须是PDB格式 (.pdb)")
        
        # 读取文件内容
        pdb_content = await file.read()
        pdb_content = pdb_content.decode('utf-8')
        
        # 验证口袋信息格式（如果提供）
        if pocket_info and pocket_info.strip():
            try:
                parse_pocket_info(pocket_info.strip())
            except ValueError as e:
                raise HTTPException(status_code=400, detail=str(e))
        
        # 处理PDB文件
        smiles, temp_file_path, message = process_pdb_to_smiles(
            pdb_content, 
            pocket_info.strip() if pocket_info and pocket_info.strip() else None
        )
        
        return PDBResponse(
            success=True,
            smiles=smiles,
            pocket_info=pocket_info.strip() if pocket_info and pocket_info.strip() else None,
            message=message,
            temp_file=temp_file_path
        )
        
    except HTTPException:
        # 重新抛出HTTP异常
        raise
    except Exception as e:
        # 清理临时文件
        if temp_file_path and os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
        
        raise HTTPException(
            status_code=500, 
            detail=f"处理PDB文件时发生错误: {str(e)}"
        )

@pdbApi.get('/health', summary='健康检查')
def health_check():
    """
    API健康检查接口
    """
    try:
        # 检查RDKit是否正常工作
        test_mol = Chem.MolFromSmiles('CCO')
        if test_mol is None:
            raise Exception("RDKit测试失败")
        
        return {
            "status": "healthy",
            "message": "PDB转SMILES服务运行正常",
            "rdkit_version": rdkit.__version__
        }
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"服务健康检查失败: {str(e)}"
        )

@pdbApi.post('/validate-pocket', summary='验证口袋信息格式')
def validate_pocket_info(pocket_info: str = Form(..., description="口袋信息字符串")):
    """
    验证口袋信息格式是否正确
    
    参数:
    - **pocket_info**: 口袋信息字符串，格式为 "A:1E8:701"
    
    返回:
    - 解析结果和验证状态
    """
    try:
        chain, ligand_name, residue_number = parse_pocket_info(pocket_info)
        return {
            "valid": True,
            "chain": chain,
            "ligand_name": ligand_name,
            "residue_number": residue_number,
            "message": "口袋信息格式正确"
        }
    except ValueError as e:
        return {
            "valid": False,
            "error": str(e),
            "message": "口袋信息格式错误"
        }