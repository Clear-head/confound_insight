from django.db import models


class Product(models.Model):
    """의약품 제품"""
    product_name = models.CharField(max_length=255, db_index=True)
    permit_number = models.CharField(max_length=50, unique=True)
    manufacturer = models.CharField(max_length=255, blank=True)
    is_combination = models.BooleanField(default=False)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = 'products'
        constraints = [
            models.UniqueConstraint(
                fields=['permit_number'],
                name='unique_product_permit_number'
            )
        ]

    def __str__(self):
        return self.product_name


class ProductIngredient(models.Model):
    """제품-성분 매핑"""

    # Choices
    STATUS_PENDING = 'PENDING'
    STATUS_SUCCESS = 'SUCCESS'
    STATUS_FAILED = 'FAILED'
    STATUS_MANUAL = 'MANUAL'

    STATUS_CHOICES = [
        (STATUS_PENDING, '대기중'),
        (STATUS_SUCCESS, '성공'),
        (STATUS_FAILED, '실패'),
        (STATUS_MANUAL, '수동'),
    ]

    product = models.ForeignKey(Product, on_delete=models.CASCADE, related_name='ingredients')
    compound = models.ForeignKey('compounds.Compound', on_delete=models.SET_NULL,
                                 null=True, blank=True, related_name='products')

    #   일단 영어로 긁어오니깐 패스
    # raw_ingredient_name = models.CharField(max_length=255)

    english_name = models.CharField(max_length=255, blank=True,
                                    help_text="MFDS MAIN_INGR_ENG - PubChem 조회용")
    is_main_active = models.BooleanField(default=True)

    normalization_status = models.CharField(max_length=20, choices=STATUS_CHOICES,
                                            default=STATUS_PENDING)
    normalization_error = models.TextField(blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = 'product_ingredients'
        constraints = [
            models.UniqueConstraint(
                fields=['product', 'english_name'],
                name='unique_product_ingredient'
            )
        ]

    def __str__(self):
        return f"{self.product.product_name}"