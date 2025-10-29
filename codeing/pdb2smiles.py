from rdkit import Chem
import re

def parse_pocket_info(pocket_str):
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

def extract_pocket_from_pdb(pdb_file, pocket_info=None):
    """
    从PDB文件中提取指定口袋区域的分子
    
    参数:
    pdb_file: PDB文件路径
    pocket_info: 口袋信息字符串，格式为 "A:1E8:701"
    
    返回:
    提取的分子对象
    """
    if pocket_info:
        # 解析口袋信息
        chain, ligand_name, residue_number = parse_pocket_info(pocket_info)
        print(f"正在提取口袋信息: 链 {chain}, 配体 {ligand_name}, 残基号 {residue_number}")
        
        # 读取PDB文件内容
        with open(pdb_file, 'r') as f:
            pdb_content = f.read()
        
        # 过滤出指定口袋区域的HETATM记录
        filtered_lines = []
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
            elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK', 'CRYST1', 'ORIGX', 'SCALE', 'MASTER', 'END')):
                # 保留PDB文件头信息
                filtered_lines.append(line)
        
        if not any(line.startswith('HETATM') for line in filtered_lines):
            print(f"警告: 未找到匹配的口袋信息 {pocket_info}")
            print("将使用完整的PDB文件进行处理")
            return Chem.rdmolfiles.MolFromPDBFile(pdb_file)
        
        # 创建临时PDB内容
        filtered_pdb_content = '\n'.join(filtered_lines)
        
        # 将过滤后的内容写入临时文件
        temp_pdb_file = pdb_file.replace('.pdb', '_pocket.pdb')
        with open(temp_pdb_file, 'w') as f:
            f.write(filtered_pdb_content)
        
        print(f"已创建口袋区域PDB文件: {temp_pdb_file}")
        
        # 从过滤后的PDB文件创建分子对象
        mol = Chem.rdmolfiles.MolFromPDBFile(temp_pdb_file)
        return mol
    else:
        # 如果没有指定口袋信息，使用完整的PDB文件
        return Chem.rdmolfiles.MolFromPDBFile(pdb_file)

def pdb_to_smiles(pdb_file, pocket_info=None):
    """
    将PDB文件转换为SMILES字符串
    
    参数:
    pdb_file: PDB文件路径
    pocket_info: 可选的口袋信息，格式为 "A:1E8:701"
    
    返回:
    SMILES字符串
    """
    try:
        # 提取分子
        mol = extract_pocket_from_pdb(pdb_file, pocket_info)
        
        if mol is None:
            print("错误: 无法从PDB文件创建分子对象")
            return None
        
        # 转换为SMILES
        smiles = Chem.MolToSmiles(mol)
        return smiles
        
    except Exception as e:
        print(f"处理过程中发生错误: {e}")
        return None

# 主程序
if __name__ == "__main__":
    pdb_file = '5p9i.pdb'
    
    # 示例1: 不使用口袋过滤
    print("=== 完整PDB文件处理 ===")
    smiles_full = pdb_to_smiles(pdb_file)
    if smiles_full:
        print(f"完整分子SMILES: {smiles_full}")
    
    # 示例2: 使用口袋信息过滤
    print("\n=== 使用口袋信息过滤 ===")
    pocket_info = "A:1E8:701"  # 链A, 配体1E8, 残基号701
    smiles_pocket = pdb_to_smiles(pdb_file, pocket_info)
    if smiles_pocket:
        print(f"口袋区域SMILES: {smiles_pocket}")
    
    print("\n=== 处理完成 ===")
