from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone
import os

# 预测结果模型 - 匹配现有数据库结构
class PredictionResult(models.Model):
    PREDICTION_TYPES = [
        ('autism_risk', '自闭症风险预测'),
    ]
    
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='用户')
    sample_id = models.IntegerField(blank=True, null=True, verbose_name='样本ID')  # 对应数据库中的 sample_id
    prediction_type = models.CharField(max_length=20, choices=PREDICTION_TYPES, default='autism_risk', verbose_name='预测类型')
    
    # 使用现有数据库字段
    confidence_score = models.FloatField(default=0.0, verbose_name='置信度分数')  # 对应数据库中的 confidence_score
    model_version = models.CharField(max_length=50, default='v1.0', verbose_name='模型版本')
    input_data = models.TextField(default='', verbose_name='输入数据')
    prediction_result = models.TextField(default='{}', verbose_name='预测结果')  # 存储JSON字符串
    
    # 时间信息
    created_time = models.DateTimeField(default=timezone.now, verbose_name='创建时间')
    
    class Meta:
        verbose_name = '预测结果'
        verbose_name_plural = '预测结果'
        ordering = ['-created_time']
        # 添加索引优化查询性能
        indexes = [
            models.Index(fields=['user', 'sample_id']),
            models.Index(fields=['created_time']),
        ]
    
    def __str__(self):
        return f"{self.user.username} - 样本{self.sample_id or '未知'} - {self.created_time}"
    
    def get_confidence(self):
        """获取置信度"""
        return self.confidence_score
    
    def get_prediction_score(self):
        """获取预测分数"""
        return self.confidence_score
    
    def get_detailed_results(self):
        """解析预测结果JSON"""
        import json
        try:
            return json.loads(self.prediction_result)
        except (json.JSONDecodeError, TypeError):
            return {}
    
    def get_risk_level(self):
        """从预测结果中获取风险等级"""
        results = self.get_detailed_results()
        return results.get('risk_level', '未知')
    
    def get_sample_name(self):
        """从预测结果中获取样本名称"""
        results = self.get_detailed_results()
        return results.get('sample_name', f'样本{self.sample_id}')

# 用户操作日志模型 - 简化版，去掉不必要的关联
class UserActivityLog(models.Model):
    ACTION_TYPES = [
        ('login', '登录'),
        ('logout', '登出'),
        ('upload', '上传文件'),
        ('analysis', '开始分析'),
        ('prediction_start', '开始预测'),
        ('view_result', '查看结果'),
        ('error', '错误'),
        ('delete', '删除数据'),
    ]
    
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='用户')
    action_type = models.CharField(max_length=20, choices=ACTION_TYPES, verbose_name='操作类型')
    description = models.CharField(max_length=500, verbose_name='操作描述')
    ip_address = models.GenericIPAddressField(verbose_name='IP地址')
    user_agent = models.TextField(verbose_name='用户代理')
    timestamp = models.DateTimeField(default=timezone.now, verbose_name='时间戳')
    
    class Meta:
        verbose_name = '用户操作日志'
        verbose_name_plural = '用户操作日志'
        ordering = ['-timestamp']
        # 添加索引优化查询性能
        indexes = [
            models.Index(fields=['user', 'action_type']),
            models.Index(fields=['timestamp']),
        ]
    
    def __str__(self):
        return f"{self.user.username} - {self.get_action_type_display()} - {self.timestamp}"
