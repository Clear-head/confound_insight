from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from django.db.backends.utils import logger

def generate_morgan_binary(smiles: str, radius: int = 2, n_bits: int = 2048) -> bytes | None:
    """
    SMILES -> Morgan Fingerprint -> Binary(bytes) 변환
    """
    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"RDKit: 유효하지 않은 SMILES입니다. ({smiles})")
            return None

        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol,
            radius=radius,
            nBits=n_bits
        )

        return fp.ToBinary()

    except Exception as e:
        logger.error(f"Fingerprint 생성 중 오류 발생: {e}")
        return None

def get_similarity(bin_fp1: bytes, bin_fp2: bytes) -> float:
    """
    저장된 두 이진 지문 간의 Tanimoto 유사도 계산 (0.0 ~ 1.0)
    """
    from rdkit.DataStructs import cDataStructs

    fp1 = cDataStructs.CreateFromBinaryText(bin_fp1.decode('latin-1'))
    fp2 = cDataStructs.CreateFromBinaryText(bin_fp2.decode('latin-1'))

    return DataStructs.TanimotoSimilarity(fp1, fp2)