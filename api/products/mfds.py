from httpx import AsyncClient
from django.conf import settings


async def search_product_by_name(search_name):
    url = "https://apis.data.go.kr/1471000/DrugPrdtPrmsnInfoService07/getDrugPrdtMcpnDtlInq07"

    all_items = []
    page_no = 1
    num_of_rows = 30

    async with AsyncClient() as client:
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
                data = response.json()

                body = data.get("body", {})
                items = body.get("items", [])
                total_count = body.get("totalCount", 0)

                all_items.extend(items)

                if len(all_items) >= total_count:
                    break

                page_no += 1

            except Exception as e:
                raise e

    return all_items
