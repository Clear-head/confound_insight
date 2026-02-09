from django.db import models


class Compound(models.Model):
    """화합물 정보"""
    standard_name = models.CharField(max_length=255, unique=True, db_index=True)
    cid = models.BigIntegerField(unique=True, null=True, blank=True)

    # 구조 정보
    smiles = models.TextField(blank=True)
    inchi = models.TextField(blank=True)
    inchi_key = models.CharField(max_length=27, blank=True)

    # 물성 정보
    molecular_formula = models.CharField(max_length=100, blank=True)
    molecular_weight = models.DecimalField(max_digits=12, decimal_places=4,
                                           null=True, blank=True)
    iupac_name = models.TextField(blank=True)

    # 지문 (RDKit)
    fingerprint_morgan = models.BinaryField(null=True, blank=True)

    # 품질 관리
    is_valid = models.BooleanField(default=True)
    validation_error = models.TextField(blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    pubchem_last_fetched = models.DateTimeField(null=True, blank=True)

    class Meta:
        db_table = 'compounds'
        constraints = [
            models.UniqueConstraint(
                fields=['standard_name'],
                name='unique_compound_standard_name'
            ),
            models.UniqueConstraint(
                fields=['cid'],
                name='unique_compound_cid',
                condition=models.Q(cid__isnull=False)
            )
        ]

    def __str__(self):
        return f"{self.standard_name} (CID:{self.cid or 'N/A'})"

    @property
    def has_structure(self):
        """구조 데이터 존재 여부"""
        return bool(self.smiles and self.fingerprint_morgan)