#!/usr/bin/env python3
import re
from collections import defaultdict
from datetime import datetime
import gzip
import pysam

hg19_contigs = [
    "##contig=<ID=chrM,length=16571>",
    "##contig=<ID=chr1,length=249250621>",
    "##contig=<ID=chr2,length=243199373>",
    "##contig=<ID=chr3,length=198022430>",
    "##contig=<ID=chr4,length=191154276>",
    "##contig=<ID=chr5,length=180915260>",
    "##contig=<ID=chr6,length=171115067>",
    "##contig=<ID=chr7,length=159138663>",
    "##contig=<ID=chr8,length=146364022>",
    "##contig=<ID=chr9,length=141213431>",
    "##contig=<ID=chr10,length=135534747>",
    "##contig=<ID=chr11,length=135006516>",
    "##contig=<ID=chr12,length=133851895>",
    "##contig=<ID=chr13,length=115169878>",
    "##contig=<ID=chr14,length=107349540>",
    "##contig=<ID=chr15,length=102531392>",
    "##contig=<ID=chr16,length=90354753>",
    "##contig=<ID=chr17,length=81195210>",
    "##contig=<ID=chr18,length=78077248>",
    "##contig=<ID=chr19,length=59128983>",
    "##contig=<ID=chr20,length=63025520>",
    "##contig=<ID=chr21,length=48129895>",
    "##contig=<ID=chr22,length=51304566>",
    "##contig=<ID=chrX,length=155270560>",
    "##contig=<ID=chrY,length=59373566>",
]


class VCFwriter:
    """
    VCF writer class
    """

    def __init__(self, vcf_out, template=None, compress=False, sample_name=None):
        self._vcf_out = vcf_out.replace(".gz", "")
        self._template = template
        self._compress = compress
        self._sample_name = sample_name

        if self._compress is True:
            self._vcf_gz_out = ("{}.gz").format(self._vcf_out)

        self._fp = open(vcf_out, "w")

        header_list = []
        if self._template is not None:
            if self.is_gz(self._template):
                f = gzip.open(self._template, "rt")
            else:
                f = open(self._template, "r")
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("##fileformat=VCF"):
                    header_list.append(line)
                    date_str = datetime.now().strftime("%d%m%y")
                    file_date = ("##filedate={}").format(date_str)
                    header_list.append(file_date)
                    continue
                if line.startswith("#"):
                    header_list.append(line)
                if line.startswith("#CHROM"):
                    tmp = line.split("\t")
                    self._sample_name = tmp[-1]
            f.close()
        else:
            header_list.append("##fileformat=VCFv4.2")
            date_time = ("##filedate={}").format(datetime.now().strftime("%d%m%y"))
            header_list.append(date_time)
            header_title = (
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}"
            ).format(self._sample_name)
            for contig in hg19_contigs:
                header_list.append(contig)
            header_list.append(header_title)

        if not self._sample_name:
            msg = "Missing sample name SN from the template vcf header"
            raise ValueError(msg)
        self._header = "\n".join(header_list)

    @staticmethod
    def is_gz(vcf) -> bool:
        with open(vcf, "rb") as test_f:
            return test_f.read(2) == b"\x1f\x8b"

    @property
    def header(self) -> str:
        """ """
        return self._header

    def add_info_to_header(self, info_dict):
        """ """
        # parse header
        keys = ["ID", "Number", "Type", "Description"]
        for key in keys:
            try:
                info_dict[key]
            except KeyError:
                msg = ("Missing {} key for INFO field").format(key)
                raise KeyError(msg)

        if not re.search('^"*"$', info_dict["Description"]):
            info_dict["Description"] = ('"{}"').format(info_dict["Description"])

        info_list = list()
        header_list = self._header.split("\n")
        empty_info = True
        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("##INFO"):
                empty_info = False
                line = line.replace("##INFO=", "")
                field_list = line.split(",")
                entry_dict = defaultdict(dict)
                for field in field_list:
                    field = field.replace("<", "").replace(">", "")
                    tmp_list = field.split("=")
                    if len(tmp_list) < 2:
                        continue
                    key = tmp_list[0]
                    val = tmp_list[1]
                    entry_dict[key] = val
                info_list.append(entry_dict)

        for entry in info_list:
            if info_dict["ID"] not in entry["ID"]:
                info_list.append(info_dict)
                break
        if empty_info is True:
            info_list.append(info_dict)

        new_header_list = list()
        seen = False

        header_list = self._header.split("\n")
        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("##"):
                if line.startswith("##INFO"):
                    continue
                new_header_list.append(line)

        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("##INFO"):
                if seen is True:
                    continue
                if seen is False:
                    for entry in info_list:
                        new_info = (
                            "##INFO=<ID={},Number={},Type={},Description={}>".format(
                                entry["ID"],
                                entry["Number"],
                                entry["Type"],
                                entry["Description"],
                            )
                        )
                        new_header_list.append(new_info)
                    seen = True

        if empty_info is True:
            for entry in info_list:
                new_info = "##INFO=<ID={},Number={},Type={},Description={}>".format(
                    entry["ID"], entry["Number"], entry["Type"], entry["Description"]
                )
                new_header_list.append(new_info)

        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("#CHROM"):
                new_header_list.append(line)

        self._header = "\n".join(new_header_list)
        return self._header

    def write(self, record):
        """ """
        self._fp.write(record + "\n")

    def close(self):
        """ """
        self._fp.close()
        if self._compress is True:
            pysam.tabix_compress(self._vcf_out, self._vcf_gz_out, force=True)
            pysam.tabix_index(self._vcf_gz_out, preset="vcf", force=True)

    def __del__(self):
        self.close()


class Record:
    """
    Variant record representation from a vcf

    :param variant_str: VCF entry as string
    :type: str
    """

    def __init__(self, variant_str: str, header=None):
        self._variant_str = variant_str
        self._header = header
        self._variant_list = self._variant_str.split("\t")
        self._custom_info = {}

        if self._header is not None:
            header_list = self._header.split("\n")
            for line in header_list:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    tmp = line.split(",")
                    id = tmp[0].replace("##INFO=<ID=", "")
                    # Now check if there's any form of complex info such as VEP CSQ
                    tmp_fmt = line.split("Format: ")

                    custom_dict = defaultdict(dict)
                    if len(tmp_fmt) > 1:
                        tmp_fmt[0] = tmp_fmt[0].replace(" ", "")
                        tmp_fmt[1] = tmp_fmt[1].replace(" ", "")
                        data = tmp_fmt[1].replace('">', "")
                        fields = data.split("|")
                        for i, field in enumerate(fields):
                            field = field.replace(" ", "")
                            custom_dict[i] = field
                        self._custom_info[id] = custom_dict
        self._INFO = self.read_info()

    def __repr__(self):
        return str(self.__dict__)

    def as_string(self) -> str:
        """ """
        tmp = []
        tmp.append(self.CHROM)
        tmp.append(self.POS)
        tmp.append(self.ID)
        tmp.append(self.REF)
        tmp.append(self.ALT)
        tmp.append(self.QUAL)
        tmp.append(self.FILTER)

        info_list = []
        subitem_info_list = []
        for key in self._INFO:
            if key in self._custom_info:
                complex_info_list = []
                item = ""
                for idx, item in enumerate(self._INFO[key]):
                    subitem_info_list = []
                    if type(item) != str:
                        for subkey in self._INFO[key][idx]:
                            subitem_info_list.append(self._INFO[key][idx][subkey])
                        subitem_info_str = "|".join(subitem_info_list)
                        complex_info_list.append(subitem_info_str)
                    else:
                        subitem_info_list.append(item)
                if type(item) != str:
                    complex_info_str = ",".join(complex_info_list).replace(" ", "_")
                else:
                    complex_info_str = "|".join(subitem_info_list).replace(" ", "_")
                value = ("{}={}").format(key, complex_info_str)
                info_list.append(value)
            else:
                if type(self._INFO[key]) is list:
                    value = self._INFO[key][0]
                else:
                    value = self._INFO[key]
                info_data = ("{}={}").format(key, value)
                info_list.append(info_data)
        tmp.append(";".join(info_list))
        tmp.append(":".join(self.FORMAT))
        tmp.append(self.SAMPLE)
        # print(self.SAMPLE)
        return "\t".join(tmp)

    def as_dict(self) -> dict():
        """ """
        var_dict = {
            "CHROM": self.CHROM,
            "POS": self.POS,
            "END": self.END,
            "ID": self.ID,
            "REF": self.REF,
            "ALT": self.ALT,
            "QUAL": self.QUAL,
            "FILTER": self.FILTER,
            "INFO": self.INFO,
            "FORMAT": self.FORMAT,
            "SAMPLE": self.SAMPLE,
        }
        return var_dict

    @property
    def CHROM(self) -> str:
        try:
            self._variant_list[0]
        except TypeError:
            raise TypeError("missing CHROM")
        else:
            return self._variant_list[0]

    @property
    def POS(self) -> str:
        try:
            self._variant_list[1]
        except TypeError:
            raise TypeError("missing POS")
        else:
            return self._variant_list[1]

    @property
    def END(self):
        end = self.POS
        if "END" in self.INFO:
            end = self.INFO["END"][0]
        return end

    @property
    def ID(self) -> str:
        alt = self._variant_list[2]
        return alt

    @property
    def REF(self) -> str:
        ref = self._variant_list[3]
        return ref

    @property
    def ALT(self) -> str:
        alt = self._variant_list[4]
        return alt

    @property
    def QUAL(self) -> str:
        qual = self._variant_list[5]
        return qual

    @property
    def FILTER(self) -> str:
        filter = self._variant_list[6]
        return filter

    @property
    def INFO(self) -> dict:
        """ """
        return self._INFO

    def read_info(self) -> dict:
        """ """
        info_dict = defaultdict(list)
        tmp_info = self._variant_list[7].split(";")
        for field in tmp_info:
            tmp_field = field.split("=", 1)
            if len(tmp_field) < 2:
                continue
            key = tmp_field[0]
            value = tmp_field[1]
            feature_list = []
            if key in self._custom_info:
                tmp_element = value.split(",")
                for element in tmp_element:
                    tmp_custom = element.split("|")
                    idx_dict = self._custom_info[key]
                    element_dict = defaultdict(dict)
                    if len(tmp_custom) > 1:
                        for i in idx_dict:
                            # print(i)
                            custom_key = idx_dict[i]
                            custom_value = tmp_custom[i]
                            if custom_key == "Consequence":
                                c_list = []
                                conseq_list = custom_value.split("&")
                                if len(conseq_list) > 1:
                                    for v in conseq_list:
                                        c_list.append(v)
                                    custom_value = "&".join(c_list)
                                else:
                                    pass
                            if "Phastcons" in custom_key or "PhyloP" in custom_key:
                                value_list = custom_value.split("&")
                                float_value_list = list()
                                new_value = custom_value
                                for item in value_list:
                                    try:
                                        float(item)
                                    except ValueError:
                                        continue
                                    else:
                                        float_value_list.append(float(item))
                                if len(float_value_list) > 0:
                                    new_value = sum(float_value_list) / len(
                                        float_value_list
                                    )
                                custom_value = str(new_value)
                            if custom_key == "HGVSp":
                                tmp_hgvsp = custom_value.split(":")
                                if len(tmp_hgvsp) > 1:
                                    custom_value = tmp_hgvsp[1]
                            if custom_key == "HGVSc":
                                tmp_hgvsc = custom_value.split(":")
                                if len(tmp_hgvsc) > 1:
                                    custom_value = tmp_hgvsc[1]
                            if not custom_value:
                                custom_value = "."
                            element_dict[custom_key] = custom_value
                    feature_list.append(element_dict)
                info_dict[key] = feature_list
            else:
                info_dict[key] = value
        return info_dict

    @INFO.setter
    def INFO(self, value_dict):
        self._INFO = value_dict

    @property
    def FORMAT(self) -> list():
        format = self._variant_list[8].split(":")
        return format

    @property
    def SAMPLE(self) -> str:
        if len(self._variant_list) > 10:
            sample = "\t".join(
                [self._variant_list[i] for i in range(9, len(self._variant_list), 1)]
            )
        else:
            sample = self._variant_list[9]
        return sample

    @property
    def len_ref(self) -> int:
        return len(self.REF)

    @property
    def GT(self) -> list:
        GT = []
        idx = len(self._variant_list) - 9
        for i in range(idx):
            GT_value = self._variant_list[8 + i + 1].split(":")
            GT_dict = {}
            for g in range(len(GT_value)):
                GT_dict[self.FORMAT[g]] = GT_value[g]
            GT.append(GT_dict)
        return GT

    @property
    def len_alt(self) -> int:
        return len(self.ALT)

    def is_snv(self) -> bool:
        """
        Return True if the variant is an SNV
        """
        if "SVTYPE" not in self.INFO:
            if self.len_ref == self.len_alt:
                if self.len_ref == 1 and self.len_alt == 1:
                    return True
                else:
                    return False
        else:
            False

    def is_mnv(self) -> bool:
        """
        Return True if the variant is an SNV
        """
        if "SVTYPE" not in self._record_dict["INFO"]:
            if self.len_ref == self.len_alt:
                if self.len_ref == 1 and self.len_alt == 1:
                    return True
                else:
                    return False
        else:
            False

    def is_deletion(self) -> bool:
        """
        Return True if the variant is a deletion (not SVTYPE)
        """
        if "SVTYPE" not in self.INFO:
            if self.len_ref > self.len_alt:
                return True
            else:
                return False
        else:
            False

    def is_insertion(self) -> bool:
        """
        Return True if the variant is an insertion
        """
        if "SVTYPE" not in self.INFO:
            if self.len_ref < self.len_alt:
                return True
            else:
                return False
        else:
            return False

    def is_cna(self) -> bool:
        """
        Return True if the variant is an insertion
        """
        # SVTYPE includes SVs and CNAs
        if "SVTYPE" in self.INFO:
            # CNAs
            if "CNA_GENES" in self.INFO:
                return True
            else:
                return False
        else:
            return False

    def is_sv(self) -> bool:
        """
        Return True if the variant is an insertion
        """
        # SVTYPE includes SVs and CNAs
        if "SVTYPE" in self.INFO:
            # CNAs
            if "CNA_GENES" in self.INFO:
                return False
            else:
                return True
        else:
            return False

    @property
    def depth(self):
        depth = "."
        if self.vartype == "SV" or self.vartype == "CNA":
            sr_ref, sr_alt = self.split_reads
            pe_ref, pe_alt = self.discordant_reads
            if sr_ref != ".":
                if sr_ref + sr_alt > 0:
                    depth = int(sr_alt) + int(sr_ref)
                elif pe_ref + pe_alt > 0:
                    depth = int(pe_alt) + int(pe_ref)
        else:
            for item in self.GT:
                if "DP" in item:
                    if item["DP"] != ".":
                        depth = item["DP"]
        return depth

    @property
    def read_support(self) -> int:
        read_support = 0
        if self.vartype == "SV" or self.vartype == "CNA":
            sr_ref, sr_alt = self.split_reads
            pe_ref, pe_alt = self.discordant_reads
            if sr_ref != ".":
                if sr_ref + sr_alt > 0:
                    read_support += sr_alt
                elif pe_ref + pe_alt > 0:
                    read_support += pe_alt
        else:
            for item in self.GT:
                if "AD" in item:
                    if item["AD"] != ".":
                        ad_list = item["AD"].split(",")
                        read_support = int(ad_list[1])
        return read_support

    @property
    def vaf(self):
        vaf = "."
        if self.vartype == "SV" or self.vartype == "CNA":
            sr_ref, sr_alt = self.split_reads
            pe_ref, pe_alt = self.discordant_reads
            if sr_ref != ".":
                if sr_ref + sr_alt > 0:
                    vaf = str(round((sr_alt / (sr_alt + sr_ref)), 3))
                elif pe_ref + pe_alt > 0:
                    vaf = str(round((pe_alt / (pe_alt + pe_ref)), 3))
        else:
            for item in self.GT:
                if "AD" in item:
                    if item["AD"] != ".":
                        ad_list = item["AD"].split(",")
                        read_support = int(ad_list[1])
                        vaf = str(round(read_support / int(self.depth), 3))
        return vaf

    @property
    def genotype(self) -> str:
        genotype = "."
        for item in self.GT:
            if "GT" in item:
                if item["GT"] != ".":
                    genotype = item["GT"]
            if "CN" in item:
                if item["CN"] != ".":
                    genotype = item["CN"]
        return genotype

    @property
    def vartype(self) -> str:
        """
        Return variant type (SNV, Deletion, Insertion, MNV, SV, CNA)
        """
        vartype = "."

        # SVTYPE includes SVs and CNAs
        if "SVTYPE" in self.INFO:
            # CNAs
            if "CNA_GENES" in self.INFO:
                vartype = "CNA"
            else:
                vartype = "SV"
        else:
            if self.len_ref == self.len_alt:
                if self.len_ref == 1 and self.len_alt == 1:
                    vartype = "SNV"
                else:
                    vartype = "MNV"
            if self.len_ref > self.len_alt:
                vartype = "Deletion"
            if self.len_ref < self.len_alt:
                vartype = "Insertion"
        return vartype

    @property
    def discordant_reads(self):
        ref = 0
        alt = 0
        for item in self.GT:
            if "PR" in item:
                if item["PR"] != ".":
                    (ref, alt) = item["PR"].split(",")
        return int(ref), int(alt)

    @property
    def split_reads(self):
        ref = 0
        alt = 0
        for item in self.GT:
            if "SR" in item:
                if item["SR"] != ".":
                    (ref, alt) = item["SR"].split(",")
        return int(ref), int(alt)

    @property
    def fold_change(self):
        fold_change = "."
        if "FOLD_CHANGE" in self.INFO:
            fold_change = self._record_dict["INFO"]["FOLD_CHANGE"]
        return fold_change

    @property
    def svend(self):
        svend = "."
        for key in self.INFO:
            if key == "END":
                svend = self.INFO["END"][0]
        return svend

    @property
    def svlen(self):
        sv_length = "."
        if self.svend != "":
            sv_length = str(int(self.svend) - int(self.POS))
        return sv_length

    @property
    def copy_number(self):
        cn = "."
        for item in self.GT:
            if "CN" in item:
                if item["CN"] != ".":
                    cn = item["CN"]
        return cn

    def set_info_subfield_from_list(self, key, value_list):
        if key not in self.INFO:
            msg = (" Missing field {} at INFO").format(key)
            raise KeyError(msg)
        self._INFO[key] = value_list

    def set_info_subfield_from_str(self, key, value_str):
        if key not in self.INFO:
            msg = (" Missing field {} at INFO").format(key)
            raise KeyError(msg)
        self._INFO[key] = value_str

    def add_info_field(self, key: str, value: str, custom: bool):
        if key not in self.INFO:
            self._INFO[key] = value
            # print(self._INFO)
        if custom is True:
            self._custom_info[key] = {}

    @property
    def civic_items(self):
        civic_items = {}
        if "CIVIC" in self.INFO:
            civic_items = self.INFO["CIVIC"]
        return civic_items

    @property
    def therapeutic_drugs(self):
        drug_list = []
        if self.civic_items:
            for item in self.civic_items:
                if "EV_DRUGS" in item:
                    if item["EV_DRUGS"]:
                        drug_results = item["EV_DRUGS"].split("&")
                        for drug in drug_results:
                            if drug != "." and drug not in drug_list:
                                drug_list.append(drug)
        return drug_list

    @property
    def therapeutic_drugs_str(self):
        drug_list = []

        if self.civic_items:
            for item in self.civic_items:
                if "EV_DRUGS" in item:
                    drug_results = item["EV_DRUGS"].split("&")
                    for drug in drug_results:
                        if drug != "." and drug not in drug_list:
                            drug_list.append(drug)
        if len(drug_list) > 0:
            return ",".join(drug_list)
        else:
            return "."

    @property
    def diseases(self):
        disease_list = []
        if self.civic_items:
            for item in self.civic_items:
                if "EV_DISEASE" in item:
                    disease_results = item["EV_DISEASE"].split("&")
                    for disease in disease_results:
                        disease = disease.replace("_", " ")
                        if disease != "." and disease not in disease_list:
                            disease_list.append(disease)
        return disease_list

    @property
    def diseases_str(self):
        disease_set = set()
        if self.civic_items:
            for item in self.civic_items:
                if "EV_DISEASE" in item:
                    disease_results = item["EV_DISEASE"].split("&")
                    for disease in disease_results:
                        disease = disease.replace("_", " ")
                        if disease != ".":
                            disease_set.add(disease)
        if len(disease_set) > 0:
            return ",".join(disease_set)
        else:
            return "."

    @property
    def clinical_trials(self):
        clinical_trials = []
        if self.civic_items:
            for item in self.civic_items:
                if "EV_CLINICAL_TRIALS" in item:
                    ctrials_list = item["EV_CLINICAL_TRIALS"].split("&")
                    for ctrial in ctrials_list:
                        if ctrial != ".":
                            clinical_trials.append(ctrial)
        return clinical_trials

    @property
    def is_rare(self) -> bool:
        is_rare = False
        if "CSQ" in self.INFO:
            if type(self.INFO["CSQ"]) is list:
                item = self.INFO["CSQ"][0]
            else:
                item = self.INFO["CSQ"]
            if item["MAX_AF_POPS"] != ".":
                if item["MAX_AF"] == ".":
                    is_rare = True
                else:
                    if float(item["MAX_AF"]) < 0.01:
                        is_rare = True
            else:
                if item["MAX_AF"] == ".":
                    is_rare = True
                else:
                    if float(item["MAX_AF"]) < 0.01:
                        is_rare = True
        return is_rare

    @property
    def clinical_trials_str(self):
        clinical_trials = []
        if self.civic_items:
            for item in self.civic_items:
                if "EV_CLINICAL_TRIALS" in item:
                    if item["EV_CLINICAL_TRIALS"]:
                        ctrials_list = item["EV_CLINICAL_TRIALS"].split("&")
                        for ctrial in ctrials_list:
                            if ctrial != ".":
                                clinical_trials.append(ctrial)
        if len(clinical_trials) > 0:
            return ",".join(clinical_trials)
        else:
            return "."

    @property
    def clinical_directions(self):
        clinical_directions = []
        if self.civic_items:
            for item in self.civic_items:
                if "EV_DIRECTION" in item:
                    if item["EV_DIRECTION"]:
                        cdirection_list = item["EV_DIRECTION"].split("&")
                        for direction in cdirection_list:
                            if direction != ".":
                                clinical_directions.append(direction)
        return clinical_directions

    @property
    def clinical_significance(self):
        clinical_significances = []

        if self.civic_items:
            for item in self.civic_items:
                if "EV_SIGNIFICANCE" in item:
                    if item["EV_SIGNIFICANCE"]:
                        signif_list = item["EV_SIGNIFICANCE"].split("&")
                        for signif in signif_list:
                            if signif != ".":
                                clinical_significances.append(signif)
        return clinical_significances

    def get_clinical_trials(self):
        clinical_trials = ""
        try:
            self.clinical_trials
        except TypeError:
            return clinical_trials
        else:
            return ",".join(self.clinical_trials)

    def get_diseases(self):
        diseases = ""
        try:
            self.diseases
        except TypeError:
            return diseases
        else:
            return ",".join(self.diseases)

    def get_therapeutic_drugs(self):
        drugs = ""
        try:
            self.drugs
        except ValueError:
            return drugs
        else:
            return ",".join(self.drugs)

    def get_tier_classification(self, source):
        """ """
        clinical_classification = None
        is_druggable = False
        is_supporting = False
        is_sensitive = False

        sensitive_significances = ["Sensitivity/Response", "Outcome", "Function"]

        for drug in self.therapeutic_drugs:
            if drug != ".":
                is_druggable = True
        for direction in self.clinical_directions:
            if "Supports" in direction:
                is_supporting = True
        for significance in self.clinical_significance:
            for s in sensitive_significances:
                if s in significance:
                    is_sensitive = True

        is_rare = False
        clinvar = "."
        impact = "."
        revel_score = "."

        if self.vartype == "SV":
            if "FUSION" in self.INFO:
                if type(self.INFO["FUSION"]) is list:
                    if self.INFO["FUSION"][0]["PARTNERS"] != ".":
                        clinical_classification = "Other"
                    else:
                        clinical_classification = "Rare"
            return clinical_classification

        if self.vartype == "CNA":
            clinical_classification = "Rare"
            return clinical_classification

        if "CSQ" in self.INFO:
            for item in self.INFO["CSQ"]:
                if item["MAX_AF_POPS"] != ".":
                    if float(item["MAX_AF"]) < 0.01:
                        is_rare = True
                else:
                    if item["MAX_AF"] == ".":
                        is_rare = True
                    else:
                        if float(item["MAX_AF"]) < 0.01:
                            is_rare = True
                clinvar = item["CLIN_SIG"]
                impact = item["IMPACT"]
                revel_score = item["REVEL_score"]
                break
        if is_rare is True:
            if is_supporting is True and is_sensitive is True and is_druggable is True:
                clinical_classification = "Therapeutic"
            elif (
                is_supporting is True and is_sensitive is False and is_druggable is True
            ):
                clinical_classification = "Other"
            else:
                if self.vartype == "SV":
                    # require paired-end reads and split-reads
                    if self.discordant_reads[1] != "." and self.split_reads[1] != ".":
                        if (
                            int(self.discordant_reads[1]) > 3
                            and int(self.split_reads[1]) > 3
                        ):
                            clinical_classification = "Other"
                elif self.vartype == "CNA":
                    pass
                else:
                    if "conflict" not in clinvar and "pathogenic" in clinvar:
                        clinical_classification = "Other"

                    if impact == "HIGH":
                        clinical_classification = "Other"
                    else:
                        if revel_score != ".":
                            if float(revel_score) > 0.75:
                                clinical_classification = "Other"
                            else:
                                # Tier 3: rare variants that do not have enough evidence
                                clinical_classification = "Rare"
                        else:
                            clinical_classification = "Rare"
        else:
            clinical_classification = "Common"

        return clinical_classification


class VCFreader:
    def __init__(self, vcf_in: str):

        self._input_vcf = vcf_in
        tmp_header = list()

        # Open gzipped file when applicable
        if self.is_gz(self._input_vcf):
            f = gzip.open(self._input_vcf, "rt")
        else:
            f = open(self._input_vcf, "r")

        # Fetch header
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                tmp_header.append(line)
            else:
                break
        f.close()
        self._header = "\n".join(tmp_header)

    def fetch(self, chr: str, pos: int = None, end: int = None) -> Record:
        """
        Fetch variants from a chr/region
        """
        tb_vcf = pysam.TabixFile(self._input_vcf)
        for rec in tb_vcf.fetch(chr, pos, end):
            record = Record(rec, self._header)
            yield record

    @staticmethod
    def is_gz(vcf: str) -> bool:
        """ """
        with open(vcf, "rb") as test_f:
            return test_f.read(2) == b"\x1f\x8b"

    @property
    def header(self) -> str:
        """ """
        return self._header

    def add_info_to_header(self, info_dict: dict) -> str:
        """ """
        # parse header
        keys = ["ID", "Number", "Type", "Description"]
        for key in keys:
            if key not in info_dict:
                msg = "Missing {} key for INFO field".format(key)
                raise KeyError(msg)

        if not re.search('^"*"$', info_dict["Description"]):
            info_dict["Description"] = ('"{}"').format(info_dict["Description"])

        info_list = list()
        header_list = self._header.split("\n")
        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("##INFO"):
                line = line.replace("##INFO=", "")
                field_list = line.split(",")
                entry_dict = defaultdict(dict)
                for field in field_list:
                    field = field.replace("<", "").replace(">", "")
                    tmp_list = field.split("=")
                    if len(tmp_list) < 2:
                        continue
                    key = tmp_list[0]
                    val = tmp_list[1]
                    entry_dict[key] = val
                info_list.append(entry_dict)

        for entry in info_list:
            if info_dict["ID"] not in entry["ID"]:
                info_list.append(info_dict)
                break

        new_header_list = list()
        seen = False

        header_list = self._header.split("\n")
        for line in header_list:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if line.startswith("##INFO"):
                    if seen is True:
                        continue
                    for entry in info_list:
                        new_info = (
                            "##INFO=<ID={},Number={},Type={},Description={}>"
                        ).format(
                            entry["ID"],
                            entry["Number"],
                            entry["Type"],
                            entry["Description"],
                        )
                        new_header_list.append(new_info)
                    seen = True
                else:
                    new_header_list.append(line)

        self._header = "\n".join(new_header_list)
        return self._header

    def __iter__(self):

        if self.is_gz(self._input_vcf):
            f = gzip.open(self._input_vcf, "rt")
        else:
            f = open(self._input_vcf, "r")
        for rec in f:
            rec = rec.rstrip("\n")
            if rec.startswith("#"):
                continue
            record = Record(rec, self._header)
            yield record
        f.close()
