cat > /root/diagnostico.sh << 'EOF'
#!/bin/bash
echo "=== DC ==="
echo ""
echo "=== 1. quality_metrics.csv ==="
head -5 lrd_workspace/analysis/quality_metrics.csv
echo "..."
wc -l lrd_workspace/analysis/quality_metrics.csv
echo ""
echo "=== 2. lrd_candidates_z4_v2.csv ==="
head -5 lrd_workspace/catalogs/lrd_candidates_z4_v2.csv
echo "..."
wc -l lrd_workspace/catalogs/lrd_candidates_z4_v2.csv
echo ""
echo "=== 3. offset_test.csv ==="
cat lrd_workspace/analysis/offset_test.csv
echo ""
echo "=== 4. fit_line_results.csv ==="
cat lrd_workspace/analysis/fit_line_results.csv
echo ""
echo "=== 5. model_comparison.csv ==="
cat lrd_workspace/analysis/model_comparison.csv
echo ""
echo "=== 6. forecast_conservative.csv ==="
cat lrd_workspace/analysis/forecast_conservative.csv
echo ""
echo "=== 7. stack_line_metrics.csv ==="
cat lrd_workspace/analysis/stack_line_metrics.csv
echo ""
echo "=== 8. lrd_classification.csv ==="
head -5 lrd_workspace/analysis/lrd_classification.csv
wc -l lrd_workspace/analysis/lrd_classification.csv
echo ""
echo "=== 9. Tier counts ==="
wc -l lrd_workspace/analysis/tier_A_files.txt
wc -l lrd_workspace/analysis/tier_B_files.txt
wc -l lrd_workspace/analysis/tier_C_files.txt
echo ""
echo "=== 10. Field counts (LRDs) ==="
echo "GOODS-S: $(ls lrd_workspace/spectra_lrd_candidates/ | grep -c 'goods-s')"
echo "GOODS-N: $(ls lrd_workspace/spectra_lrd_candidates/ | grep -c 'goods-n')"
echo ""
echo "=== 11. Field counts (z>5) ==="
echo "GOODS-S: $(ls lrd_workspace/spectra/ | grep -c 'goods-s')"
echo "GOODS-N: $(ls lrd_workspace/spectra/ | grep -c 'goods-n')"
echo ""
echo "=== 12. Field counts (todos) ==="
echo "GOODS-S: $(ls todos_los_fits/ | grep -c 'goods-s')"
echo "GOODS-N: $(ls todos_los_fits/ | grep -c 'goods-n')"
echo ""
echo "=== 13. master_features_v4.csv ==="
head -2 lrd_workspace/analysis/master_features_v4.csv
echo "..."
wc -l lrd_workspace/analysis/master_features_v4.csv
echo ""
echo "=== FIN ==="
EOF
chmod +x /root/diagnostico.sh
/root/diagnostico.sh > /root/diagnostico_output.txt
cat /root/diagnostico_output.txt