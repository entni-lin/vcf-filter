# VCF Filtering Tool

VCF Filtering Tool is a Python command-line tool designed to filter Mutect2-generated VCF files using user-defined
criteria provided in a JSON file.

## Requirements
- **Python Version:** tested with Python 3.12.4
- **External Modules:**  
  - [pysam](https://github.com/pysam-developers/pysam) — for efficient VCF parsing and streaming. To install the module, run:
    ```bash
    pip install pysam
    ```

## Usage
Run the tool from the project's root directory using the following command:
```bash
python3 scripts/vcf_filter.py --vcf data/input.vcf --criteria config/criteria.json --output data/output.vcf
```
- --vcf: Specifies the input VCF file in Mutect2 format (e.g., data/input.vcf).
- --criteria: Points to the JSON file defining the filtering criteria. For example:
```json
{
    "TLOD": ">=10",
    "DP": ">=20"
}
```
- **Note:** If your raw VCF data does not set any values for the FILTER field (e.g., artifact), it is advisable to remove any FILTER-based comparisons from your JSON criteria.
- --output: The output file where the filtered VCF will be written (*e.g.,* data/output.vcf).

To inspect the output, users can use the following commands:

- Preview the first 20 lines in an organized table (with headers):
```bash
head -n 20 data/output.vcf | grep -v '^##' | column -t -s $'\t'
```
- Display only the FILTER column (the 7th column):
```bash
grep -v '^##' data/output.vcf | cut -f7 | head -n 20
```

## Testing
To ensure that the tool works as expected, consider the following testing steps:

1. **Prepare Test Data:** In the `data/` folder, include one or more sample VCF files (*e.g.,* `data/SAMPLE_mutect2_raw.vcf`) along with a corresponding JSON criteria file (*e.g.,* `data/filter_criteria.json`).
2. **Run the tool:** Install the required module `pysam` (if not already installed) and run the command-line tool. This will process each variant in users' input VCF file and output a new file where variants meeting the criteria have their FILTER column updated to `PASS`.
3. **Manual Validation:** After running the tool, use the commands above to check whether variants meeting the criteria have their FILTER columns updated to `PASS`.
4. **Compare Expected Output (Optional):** If an expected output file is available (*e.g.,* data/expected_output.vcf), use the diff command to compare it with the generated output:
```
diff data/output.vcf data/expected_output.vcf
```
 - No differences indicate that the tool is functioning properly.

## Design Decisions

- **VCF Parsing with pysam for streaming processing:**  
  We chose pysam for VCF parsing because it supports streaming processing, allowing the tool to handle very large VCF files without loading the entire file into memory. This not only improves efficiency but also reduces memory usage.

- **Handling Multi-Allelic Sites:**  
  For variant records that may contain multiple alternate alleles (*i.e.,* multi-allelic sites) in the INFO field, the tool is designed to iterate over lists/tuples. In the filtering logic, if any one of the values in the multi-allelic field meets the specified criteria, the record is considered to pass the filtering stage. This design strikes a balance between complexity and functionality while accommodating common VCF scenarios.

- **Criteria Parsing and Flexibility:**  
  Filtering criteria are provided via a JSON file. The tool parses each criterion into an operator and a threshold value. This design decision allows for easy extension of logic; new conditional operators or custom filtering rules can be added later with minimal code changes.

- **Modular Code Structure:**  
  The module is structured into functions (*e.g.,* one function for argument parsing, one for criteria evaluation, and another for VCF processing). This modular approach enhances maintainability and simplifies future debugging or extension efforts.

## Design Concept and Workflow

**Core Algorithm**  

The tool follows a structured workflow:

1. **Argument Parsing:**  
It uses the `argparse` module to accept input parameters such as the VCF file, JSON criteria file, and output file.
2. **Criteria Loading and Parsing:**  
The criteria in the JSON file are read and parsed into operator-threshold pairs. This allows the tool to dynamically compare various fields (like TLOD and DP) against the specified conditions.
3. **VCF Processing Using Streaming:**  
Utilizing the `pysam.VariantFile` method, the VCF file is streamed record by record. This streaming approach ensures that even large VCF files are processed efficiently without loading the entire file into memory.
4. **Record Evaluation:**  
Each variant (record) is evaluated based on the criteria. For each record, its INFO (and, if applicable, FILTER) fields are checked. If every condition is satisfied, the FILTER field is cleared and updated to `PASS`.
5. **Output Generation:**  
Processed variant records are then written to a new output VCF file. This step is handled in a streaming fashion to maintain low resource usage.

---

### How to utilize the tool

Users can run the tool from the terminal. Ensure that Python 3.12.4 or greater is installed along with the external dependency `pysam`. 

Then, execute the program with the appropriate command-line arguments:
```bash
python3 scripts/vcf_filter.py --vcf data/input.vcf --criteria config/criteria.json --output data/output.vcf
```
This command will process the input VCF file according to the criteria defined in the JSON file and generate an output VCF file with variants that pass the criteria marked as `PASS`.

## Maintenance and Future Extensions
- Modular Code Structure: The module is organized into distinct functions. This modular design makes it easier to isolate and fix bugs, update functions independently, or add new features.
- Extending Functionality: Future improvements can include:
  1. Adding support for more complex logical operations between criteria.
  2. Enhancing the filtering criteria to handle additional comparison operators or more intricate filtering logic.
  3. Refining multi-allelic site handling for more complex datasets.
- Version Control: With this module maintained under a version control system (*e.g.,* Git), managing changes, collaborating with other developers, and rolling back to previous versions becomes straightforward.

## Usage Limitations
- **Limited Comparison Operators:** The tool currently supports basic numeric and string comparisons. More advanced logic (e.g., complex Boolean conditions) is not implemented.
- **Handling of Multi-allelic Sites:** When dealing with fields that may contain multiple values (e.g., due to multi-allelic sites), the tool passes a record if at least one of the values meets the filtering criterion. This approach may need refining for datasets that require more granular handling.
- **Filtering the FILTER Field:** The current implementation supports only equality checks for the FILTER field. Users whose data have different conventions for this field should adjust the JSON filter criteria accordingly.
