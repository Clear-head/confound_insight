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

    raw_ingredient_name = models.CharField(max_length=255)
    content = models.CharField(max_length=100, blank=True)
    unit = models.CharField(max_length=20, blank=True)
    is_main_active = models.BooleanField(default=True)

    normalization_status = models.CharField(max_length=20, choices=STATUS_CHOICES,
                                            default=STATUS_PENDING)
    normalization_error = models.TextField(blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = 'product_ingredients'
        unique_together = [['product', 'raw_ingredient_name']]

    def __str__(self):
        return f"{self.product.product_name} - {self.raw_ingredient_name}"