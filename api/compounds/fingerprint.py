from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from django.db.backends.utils import logger

async def generate_morgan_binary(smiles: str, radius: int = 2, n_bits: int = 2048) -> bytes | None:
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