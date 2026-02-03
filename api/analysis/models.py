from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator


class SimilarityAnalysis(models.Model):
    """화합물 간 Tanimoto 유사도 분석 결과"""
    target_compound = models.ForeignKey(
        'compounds.Compound',
        on_delete=models.CASCADE,
        related_name='similarities_as_target'
    )
    similar_compound = models.ForeignKey(
        'compounds.Compound',
        on_delete=models.CASCADE,
        related_name='similarities_as_comparison'
    )

    similarity_score = models.FloatField(
        validators=[MinValueValidator(0.0), MaxValueValidator(1.0)]
    )

    analysis_date = models.DateTimeField(auto_now_add=True)
    is_current = models.BooleanField(default=True)

    class Meta:
        db_table = 'compound_similarities'
        unique_together = [['target_compound', 'similar_compound']]
        constraints = [
            models.CheckConstraint(
                condition=~models.Q(target_compound=models.F('similar_compound')),
                name='no_self_similarity'
            )
        ]

    def __str__(self):
        return f"{self.target_compound.standard_name} ↔ {self.similar_compound.standard_name}: {self.similarity_score:.3f}"