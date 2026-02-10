from httpx import AsyncClient
from django.conf import settings
from django.db.backends.utils import logger

async def search_product_by_name(search_name: str) -> dict:
    url = "https://apis.data.go.kr/1471000/DrugPrdtPrmsnInfoService07/getDrugPrdtMcpnDtlInq07"

    result = {}
    page_no = 1
    num_of_rows = 100

    async with AsyncClient(timeout=10.0) as client:
        while True:
            params = {
                "serviceKey": settings.MFDS_API_KEY,
                "type": "json",
                "Prduct": search_name,
                "pageNo": page_no,
                "numOfRows": num_of_rows
            }

            try:
                response = await client.get(url, params=params)
                response.raise_for_status()
                res_data = response.json()

                items = res_data.get("body", {}).get("items", [])
                total_count = res_data.get("body", {}).get("totalCount", 0)


                if not items:
                    break

                for item in items:
                    item_seq = item.get("ITEM_SEQ")
                    if result.get(item_seq, None) is None:
                        result[item_seq] = {
                            "product_name": item.get("PRDUCT"),
                            "manufacturer": item.get("ENTRPS"),
                            "permit_number": item.get("MTRAL_CODE"),
                            "ingredients_en": list(map(str, item.get("MAIN_INGR_ENG").split("/"))),
                            "is_combination": False if len(item.get("MAIN_INGR_ENG").split("/")) == 1 else True,
                        }
                    else:
                        pass

                if len(items) < num_of_rows or (page_no * num_of_rows) >= total_count:
                    break
                page_no += 1

            except Exception as e:
                logger.error(f"식약처 API 데이터 처리 중 오류: {e}")
                break

    return result