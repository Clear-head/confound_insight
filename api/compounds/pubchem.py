import asyncio

from django.db.backends.utils import logger
from httpx import AsyncClient, HTTPError
from pydantic import BaseModel

base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
limit_semaphore = asyncio.Semaphore(4)
class PubChemInformation(BaseModel):
    CID: int
    MolecularFormula: str
    MolecularWeight: float
    ConnectivitySMILES: str
    InChI: str
    InChIKey: str
    IUPACName: str

async def get_cid_by_name(name: str)->int|None:

    url = base_url + f"compound/name/{name}/cids/JSON"

    async with limit_semaphore:
        async with AsyncClient() as client:
            try:

                response = await client.get(url)
                response.raise_for_status()
                res = response.json().get("IdentifierList", {}).get("CID", [])

                if not res:
                    logger.warning(f"CID를 찾을 수 없습니다: {name}")
                    return None

                if len(res) > 1:
                    logger.info(f"여러 CID가 발견되었습니다 ({len(res)}개). 첫 번째 CID 사용: {res[0]} (검색어: {name})")

                return res[0]

            except HTTPError as http_err:
                logger.error(f"HTTP error occurred: {http_err}")
                return None

            except Exception as err:
                logger.error(f"Other error occurred: {err}")
                return None

async def get_compound_by_cid(cid: int)->PubChemInformation|None:
    option = "CanonicalSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName"
    url = base_url + f"compound/cid/{cid}/property/{option}/JSON"

    async with limit_semaphore:
        async with AsyncClient() as client:
            try:
                response = await client.get(url)
                response.raise_for_status()
                res = response.json().get("PropertyTable", {}).get("Properties", [])

                if len(res) == 0:
                    raise Exception
                else:
                    return PubChemInformation(**res[0])

            except HTTPError as http_err:
                logger.error(f"HTTP error occurred: {http_err}")
                return None

            except Exception as err:
                logger.error(f"Other error occurred: {err}")
                return None