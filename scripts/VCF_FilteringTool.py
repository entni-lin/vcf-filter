"""
VCF Filtering Tool

解析來自 Mutect2 的 VCF 檔案，依據 JSON 定義的條件進行篩選，若符合所有條件則將 FILTER 欄位標記為 PASS，
否則保留原有的標記（或使用使用者自訂的標記）。

使用方式:
    python3 vcf_filter.py --vcf input.vcf --criteria criteria.json --output output.vcf
"""

import argparse # 用來解析命令列參數（command-line arguments）: 使用裡面的 class & methods
import json
import re # 正規表示式 - sequence of characters that specifies a match pattern in tex
import sys
import pysam # For VCF parsing

# 允許的比較運算子
ALLOWED_OPERATORS = {">", ">=", "<", "<=", "==", "!="}

def parse_args():
    parser = argparse.ArgumentParser(description="VCF Filtering Tool (Mutect2 format)")
    parser.add_argument("--vcf", required=True, help="輸入的 VCF 檔案 (Mutect2 format)")
    parser.add_argument("--criteria", required=True, help="定義篩選標準的 JSON 檔案 (例如 {\"TLOD\":\">=10\", \"DP\":\">=20\"})")
    parser.add_argument("--output", required=True, help="輸出的 VCF 檔案")
    return parser.parse_args()


def load_criteria(json_file):
    """
    載入 JSON 格式的條件，並解析運算子與數值
    回傳 dictionary: { field: (operator, threshold or pattern) }
    若條件字串未包含運算子，預設為相等 (==) 比較
    """
    try:
        with open(json_file, "r") as f:
            raw_criteria = json.load(f)
    except Exception as e:
        sys.exit(f"Error loading JSON file: {e}")

    criteria = {}
    # 解析每一個條件值，如 ">=10" 拆解為 (" >=", "10")
    for key, condition in raw_criteria.items():
        # 檢查 condition 是否包含明確的運算子
        match = re.match(r"^(>=|<=|==|!=|>|<)(.+)$", condition.strip()) # 若有match,回傳一個tuple e.g., ('>=', '10')
        if match:
            op, threshold = match.groups() # 拆開 運算子 和 數值
            if op not in ALLOWED_OPERATORS:
                sys.exit(f"不允許的運算子 {op} 在條件 {key}: {condition}")
            criteria[key] = (op, threshold.strip())
        else:
            # 若未指定運算子則使用 == 進行比對 (例如: FILTER == "artifact")
            criteria[key] = ("==", condition.strip())
    return criteria

def evaluate_condition(value, op, threshold):
    """
    比較 value 與 threshold 是否符合運算子（支援數字與字串比較）
    注意: 如果無法轉型為 float 則會以字串方式比較
    """
    # 先嘗試以數字方式比較
    try:
        # 同時轉換 value 與 threshold 為 float
        value_num = float(value)
        threshold_num = float(threshold)
        if op == ">":
            return value_num > threshold_num
        elif op == ">=":
            return value_num >= threshold_num
        elif op == "<":
            return value_num < threshold_num
        elif op == "<=":
            return value_num <= threshold_num
        elif op == "==":
            return value_num == threshold_num
        elif op == "!=":
            return value_num != threshold_num
    except (TypeError, ValueError):
        # 無法轉成數字時，利用字串比較
        if op == "==":
            return str(value) == threshold
        elif op == "!=":
            return str(value) != threshold
        else:
            # 對字串進行大小比較意義不大，這裡視為不符合條件
            return False

def variant_passes(record, criteria):
    """
    檢查單筆 VCF record 是否符合所有 criteria
    - 從 record.info 中獲取指定欄位 (若欄位不存在則視為不符合條件)
    - 如果條件為 FILTER，則檢查 record.filter.keys() 或 record.filter (取決於 VCF 格式)
    - 能處理數字與字串型態的條件比較
    """
    for field, (op, thresh) in criteria.items():
        # 特殊處理 FILTER 欄位: 直接取 record.filter (可能有多個標籤)
        if field.upper() == "FILTER":
            # 許多 VCF 會以一組標籤表示，多個標籤時可能存成 list 或集合
            current_filters = record.filter if hasattr(record, "filter") else record.info.get("FILTER", None)
            # 若 current_filters 是集合或 list，則判斷其中是否有任一項符合條件
            # 為簡化邏輯，這裡只檢查字串相等條件，其他運算子不常用於 FILTER 的狀況
            if op != "==":
                sys.exit("目前僅支援 FILTER 欄位的字串相等比對")
            # 轉成 list 處理
            if current_filters is None:
                return False
            if isinstance(current_filters, (list, tuple, set)):
                if thresh not in current_filters:
                    return False
            else:
                if str(current_filters) != thresh:
                    return False
        else:
            # 從 INFO 欄位取得數值或字串
            value = record.info.get(field, None)
            if value is None:
                # 若欄位缺失即視為不符合
                if field.upper() == "TLOD" or field.upper() == "DP":
                    sys.stderr.write(
                        f"Warning: Field '{field}' is missing in record at {record.chrom}:{record.pos}. Condition is not met.\n"
                    )
                return False
            # 注意: value 有時候可能為 list (例如多個等位基因的情況)
            if isinstance(value, (list, tuple)):
                # 若列表中至少有一個值符合就算通過此欄位
                if not any(evaluate_condition(val, op, thresh) for val in value):
                    return False
            else:
                if not evaluate_condition(value, op, thresh):
                    return False
    return True

def main():
    args = parse_args()
    criteria = load_criteria(args.criteria)
    try:
        # 使用 pysam 以 streaming 方式讀取 VCF 檔案，適用於大型檔案
        vcf_in = pysam.VariantFile(args.vcf, "r")
    except Exception as e:
        sys.exit(f"Error openning VCF file: {e}")

    # 建立輸出 VCF，複製 header
    try:
        vcf_out = pysam.VariantFile(args.output, "w", header=vcf_in.header)
    except Exception as e:
        sys.exit(f"Error creating output VCF file: {e}")

    # 遍歷每一筆 variant 記錄
    for record in vcf_in:
        if variant_passes(record, criteria):
            # 若 variant 符合所有條件，則將 FILTER 欄位設為 PASS（清除舊的標籤）
            record.filter.clear()  # 若 FILTER 原本有多個標籤，先清空
            record.filter.add("PASS")
        # 若不符合則保留原有的 FILTER 標籤（或可依需改成其他標籤）
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    sys.exit(0)

if __name__ == "__main__":
    main()