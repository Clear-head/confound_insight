from typing import List, Optional

from django.db.backends.utils import logger
from django.utils import timezone
from django.db import transaction

from api.compounds.models import Compound
from api.compounds.pubchem import get_cid_by_name, get_compound_by_cid
from api.compounds.fingerprint import generate_morgan_binary, get_similarity
from api.analysis.models import SimilarityAnalysis


async def get_or_create_compound(standard_name: str) -> Compound:
    """
    표준화된 성분명으로 Compound 조회 또는 생성
    """
    compound, created = Compound.objects.aget_or_create(standard_name=standard_name)
    if created:
        logger.info(f"Created compound {standard_name}")
        await enrich_compound_from_pubchem(compound)

    return compound


async def enrich_compound_from_pubchem(compound: Compound) -> bool:
    """
    PubChem API를 통해 화합물 정보를 수집하고 DB에 저장

    처리 단계:
    1. CID가 없으면 get_cid_by_name() 호출하여 조회
    2. CID로 get_compound_by_cid() 호출
    3. SMILES, 분자식, 분자량 등을 Compound 모델에 저장
    4. Fingerprint 생성 (generate_morgan_binary)
    5. pubchem_last_fetched 타임스탬프 업데이트
    """
    try:
        if not compound.cid:
            cid = await get_cid_by_name(compound.standard_name)
        else:
            cid = compound.cid
        raw_data = await get_compound_by_cid(cid)

        compound.cid = cid
        compound.molecular_formula = raw_data.MolecularFormula
        compound.molecular_weight = raw_data.MolecularWeight
        compound.smiles = raw_data.ConnectivitySMILES
        compound.inchi = raw_data.InChI
        compound.inchi_key = raw_data.InChIKey

        compound.fingerprint_morgan = generate_morgan_binary(compound.smiles)

        compound.pubchem_last_fetched = timezone.now()

        await compound.asave()

        return True

    except Exception as err:
        logger.error(f"Failed to enrich {compound.standard_name}")
        return False


async def update_compound_fingerprint(compound: Compound) -> bool:
    """
    SMILES 문자열로부터 Morgan Fingerprint를 생성하고 저장
    """
    if not compound.smiles:
        return False

    try:
        compound.fingerprint_morgan = generate_morgan_binary(compound.smiles)
        compound.pubchem_last_fetched = timezone.now()
        await compound.asave()
        return True

    except Exception as err:
        logger.error(f"Failed to enrich {compound.standard_name}")
        return False