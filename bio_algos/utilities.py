"""
Advanced Bioinformatics Toolkit
A comprehensive library for DNA/RNA sequence analysis, manipulation, and retrieval.
"""

from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
from Bio import Align, Entrez, SeqIO
import io
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import re
import os
from enum import Enum
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SequenceType(Enum):
    """Enumeration for sequence types."""
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"


class GeneticCode(Enum):
    """Enumeration for genetic code tables."""
    STANDARD = 1
    VERTEBRATE_MITOCHONDRIAL = 2
    YEAST_MITOCHONDRIAL = 3
    BACTERIAL = 11


@dataclass
class AlignmentResult:
    """Container for alignment results."""
    alignment: MultipleSeqAlignment
    score: float
    identity: float
    gaps: int


@dataclass
class SequenceAnalysis:
    """Container for sequence analysis results."""
    gc_content: float
    at_content: float
    length: int
    molecular_weight: float
    melting_temp: float
    sequence_type: SequenceType


class InvalidSequenceException(Exception):
    """Raised when an invalid sequence is provided."""
    pass


class InvalidQueryException(Exception):
    """Raised when a database query fails."""
    def __init__(self, query: str, database: str):
        self.query = query
        self.database = database
        super().__init__(f"Query '{query}' failed in database '{database}'")


class CodonTable:
    """Manages codon tables for different genetic codes."""
    
    STANDARD_DNA = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    
    STANDARD_RNA = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    
    # Start and stop codons
    PROKARYOTIC_START_CODONS = ["ATG", "GTG", "TTG"]  # DNA
    EUKARYOTIC_START_CODONS = ["ATG"]  # DNA
    STOP_CODONS_DNA = ["TAG", "TAA", "TGA"]
    STOP_CODONS_RNA = ["UAG", "UAA", "UGA"]


class SequenceTools:
    """Main class for sequence manipulation and analysis."""
    
    def __init__(self):
        self.codon_table = CodonTable()
    
    @staticmethod
    def validate_sequence(sequence: str, seq_type: SequenceType) -> bool:
        """
        Validate if a sequence contains only valid characters.
        
        Args:
            sequence: The sequence to validate
            seq_type: The type of sequence (DNA, RNA, or PROTEIN)
            
        Returns:
            bool: True if valid, False otherwise
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if seq_type == SequenceType.DNA:
            pattern = re.compile(r'^[ATCGN]+$')
        elif seq_type == SequenceType.RNA:
            pattern = re.compile(r'^[AUCGN]+$')
        elif seq_type == SequenceType.PROTEIN:
            pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY*]+$')
        else:
            return False
            
        return bool(pattern.match(sequence))
    
    def dna_to_rna(self, dna_sequence: str) -> str:
        """
        Convert DNA sequence to RNA sequence (transcription).
        
        Args:
            dna_sequence: Input DNA sequence
            
        Returns:
            str: RNA sequence
            
        Raises:
            InvalidSequenceException: If the DNA sequence is invalid
        """
        dna_sequence = dna_sequence.upper().replace('\n', '').replace(' ', '')
        
        if not self.validate_sequence(dna_sequence, SequenceType.DNA):
            raise InvalidSequenceException("Invalid DNA sequence. Only A, T, C, G, and N are allowed.")
        
        # Simple transcription: replace T with U
        return dna_sequence.replace('T', 'U')
    
    def reverse_complement(self, sequence: str, seq_type: SequenceType = SequenceType.DNA) -> str:
        """
        Get the reverse complement of a DNA or RNA sequence.
        
        Args:
            sequence: Input sequence
            seq_type: Type of sequence (DNA or RNA)
            
        Returns:
            str: Reverse complement sequence
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if seq_type == SequenceType.DNA:
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        elif seq_type == SequenceType.RNA:
            complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        else:
            raise InvalidSequenceException("Can only complement DNA or RNA sequences")
        
        if not self.validate_sequence(sequence, seq_type):
            raise InvalidSequenceException(f"Invalid {seq_type.value} sequence")
        
        # Get complement and reverse
        complement = ''.join(complement_map.get(base, 'N') for base in sequence)
        return complement[::-1]
    
    def gc_content(self, sequence: str) -> float:
        """
        Calculate GC content of a sequence.
        
        Args:
            sequence: DNA or RNA sequence
            
        Returns:
            float: GC content as a percentage (0-100)
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if len(sequence) == 0:
            return 0.0
        
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def melting_temperature(self, sequence: str, salt_conc: float = 0.05) -> float:
        """
        Calculate melting temperature using the Wallace rule for short sequences
        and salt-adjusted formula for longer sequences.
        
        Args:
            sequence: DNA sequence
            salt_conc: Salt concentration in M (default 0.05M)
            
        Returns:
            float: Melting temperature in Celsius
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if len(sequence) < 14:
            # Wallace rule: Tm = 4(G+C) + 2(A+T)
            gc = sequence.count('G') + sequence.count('C')
            at = sequence.count('A') + sequence.count('T')
            return 4 * gc + 2 * at
        else:
            # Salt-adjusted formula
            gc_content = self.gc_content(sequence) / 100
            return 81.5 + 16.6 * (0.41 + 0.025 * salt_conc) + 0.41 * gc_content - 500 / len(sequence)
    
    def translate(self, sequence: str, 
                  genetic_code: GeneticCode = GeneticCode.STANDARD,
                  to_stop: bool = True,
                  complete: bool = False) -> str:
        """
        Translate a DNA or RNA sequence to protein.
        
        Args:
            sequence: DNA or RNA sequence
            genetic_code: Which genetic code to use
            to_stop: Stop at first stop codon
            complete: Only translate complete codons
            
        Returns:
            str: Protein sequence
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        # Determine if DNA or RNA
        if 'U' in sequence:
            codon_table = self.codon_table.STANDARD_RNA
            seq_type = SequenceType.RNA
        else:
            codon_table = self.codon_table.STANDARD_DNA
            seq_type = SequenceType.DNA
        
        if not self.validate_sequence(sequence, seq_type):
            raise InvalidSequenceException(f"Invalid {seq_type.value} sequence")
        
        protein = []
        
        # Process codons
        for i in range(0, len(sequence) - 2, 3):
            if i + 3 > len(sequence) and complete:
                break
                
            codon = sequence[i:i+3]
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, 'X')  # X for unknown
                
                if amino_acid == '*' and to_stop:
                    break
                    
                protein.append(amino_acid)
        
        return ''.join(protein)
    
    def find_orfs(self, sequence: str, min_length: int = 100) -> List[Tuple[int, int, str]]:
        """
        Find all open reading frames (ORFs) in a sequence.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of tuples: (start_pos, end_pos, protein_sequence)
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        orfs = []
        
        # Check all three reading frames in both directions
        for strand, seq in [('+', sequence), ('-', self.reverse_complement(sequence))]:
            for frame in range(3):
                # Find all start codons
                for i in range(frame, len(seq) - 2, 3):
                    codon = seq[i:i+3]
                    if codon in self.codon_table.EUKARYOTIC_START_CODONS:
                        # Translate until stop codon
                        protein = self.translate(seq[i:], to_stop=True)
                        orf_length = len(protein) * 3
                        
                        if orf_length >= min_length:
                            if strand == '+':
                                orfs.append((i, i + orf_length, protein))
                            else:
                                # Adjust positions for reverse strand
                                orfs.append((len(sequence) - i - orf_length, 
                                           len(sequence) - i, protein))
        
        return orfs
    
    def analyze_sequence(self, sequence: str, seq_type: SequenceType) -> SequenceAnalysis:
        """
        Perform comprehensive sequence analysis.
        
        Args:
            sequence: Input sequence
            seq_type: Type of sequence
            
        Returns:
            SequenceAnalysis object with results
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if not self.validate_sequence(sequence, seq_type):
            raise InvalidSequenceException(f"Invalid {seq_type.value} sequence")
        
        gc_content = self.gc_content(sequence)
        at_content = 100 - gc_content
        length = len(sequence)
        
        # Calculate molecular weight (simplified)
        if seq_type == SequenceType.DNA:
            avg_mw_per_base = 330  # Average MW of DNA nucleotide
        elif seq_type == SequenceType.RNA:
            avg_mw_per_base = 340  # Average MW of RNA nucleotide
        else:
            avg_mw_per_base = 110  # Average MW of amino acid
        
        molecular_weight = length * avg_mw_per_base
        
        # Melting temperature (only for DNA)
        melting_temp = self.melting_temperature(sequence) if seq_type == SequenceType.DNA else 0.0
        
        return SequenceAnalysis(
            gc_content=gc_content,
            at_content=at_content,
            length=length,
            molecular_weight=molecular_weight,
            melting_temp=melting_temp,
            sequence_type=seq_type
        )


class SequenceAlignment:
    """Class for sequence alignment operations."""
    
    @staticmethod
    def align_sequences(sequences: List[str], 
                       ids: List[str],
                       seq_type: str = "nucleotide",
                       algorithm: str = "needleman-wunsch") -> AlignmentResult:
        """
        Align sequences using various algorithms.
        
        Args:
            sequences: List of sequences to align
            ids: List of sequence identifiers
            seq_type: Type of sequences ("nucleotide", "protein", "genome")
            algorithm: Alignment algorithm to use
            
        Returns:
            AlignmentResult object
        """
        if len(sequences) < 2:
            raise ValueError("At least two sequences are required for alignment")
        
        if len(sequences) != len(ids):
            raise ValueError("Number of sequences must match number of IDs")
        
        # Create aligner
        aligner = Align.PairwiseAligner()
        
        # Configure scoring based on sequence type
        if seq_type == 'nucleotide':
            aligner.match_score = 3
            aligner.mismatch_score = -3
            aligner.open_gap_score = -7
            aligner.extend_gap_score = -2
        elif seq_type == 'protein':
            aligner.match_score = 1
            aligner.mismatch_score = -1
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -1
        elif seq_type == 'genome':
            aligner.match_score = 1
            aligner.mismatch_score = -1
            aligner.open_gap_score = -5
            aligner.extend_gap_score = -0.5
        
        # Set algorithm
        if algorithm == "needleman-wunsch":
            aligner.mode = 'global'
        elif algorithm == "smith-waterman":
            aligner.mode = 'local'
        
        # Perform alignment
        alignments = aligner.align(sequences[0], sequences[1])
        best_alignment = next(alignments)
        
        # Create SeqRecord objects
        seq_records = []
        for i, (seq, seq_id) in enumerate(zip([str(best_alignment[0]), str(best_alignment[1])], ids[:2])):
            seq_records.append(SeqRecord(Seq(seq), id=seq_id))
        
        # Create MultipleSeqAlignment
        msa = MultipleSeqAlignment(seq_records)
        
        # Calculate statistics
        identity = SequenceAlignment._calculate_identity(str(best_alignment[0]), str(best_alignment[1]))
        gaps = str(best_alignment[0]).count('-') + str(best_alignment[1]).count('-')
        
        return AlignmentResult(
            alignment=msa,
            score=best_alignment.score,
            identity=identity,
            gaps=gaps
        )
    
    @staticmethod
    def _calculate_identity(seq1: str, seq2: str) -> float:
        """Calculate percentage identity between two aligned sequences."""
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
        return (matches / len(seq1)) * 100


class NCBITools:
    """Class for interacting with NCBI databases."""
    
    def __init__(self, email: str, api_key: Optional[str] = None, cache_dir: str = "ncbi_cache"):
        """
        Initialize NCBI tools with credentials.
        
        Args:
            email: Email address for NCBI
            api_key: Optional API key for increased rate limits
        """
        self.email = email
        self.api_key = api_key or os.environ.get('NCBI_API_KEY')
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
    
    def search_database(self, 
                       query: str,
                       database: str = "nucleotide",
                       retmax: int = 10,
                       field: Optional[str] = None) -> List[str]:
        """
        Search NCBI database and return list of IDs.
        
        Args:
            query: Search query
            database: NCBI database to search
            retmax: Maximum number of results
            field: Specific field to search
            
        Returns:
            List of record IDs
        """
        try:
            safe_query = re.sub(r"[^0-9a-zA-Z_-]", "_", query)
            field_part = f"_{field}" if field else ""
            cache_file = os.path.join(
                self.cache_dir,
                f"search_{database}_{safe_query}_{retmax}{field_part}.json",
            )
            if os.path.exists(cache_file):
                with open(cache_file, "r") as f:
                    return json.load(f)

            search_params = {
                "db": database,
                "term": query,
                "retmax": retmax,
            }
            if field:
                search_params["field"] = field

            handle = Entrez.esearch(**search_params)
            results = Entrez.read(handle)
            handle.close()

            ids = results.get("IdList", [])
            with open(cache_file, "w") as f:
                json.dump(ids, f)

            return ids

        except Exception as e:
            logger.error(f"Search failed: {e}")
            raise InvalidQueryException(query, database)
    
    def fetch_records(self, 
                     ids: Union[str, List[str]],
                     database: str = "nucleotide",
                     rettype: str = "fasta",
                     parse: bool = True) -> Union[List[SeqRecord], List[str]]:
        """
        Fetch records from NCBI database.
        
        Args:
            ids: Single ID or list of IDs
            database: NCBI database
            rettype: Return type (fasta, gb, etc.)
            parse: Whether to parse results into SeqRecord objects
            
        Returns:
            List of SeqRecord objects or raw text
        """
        if isinstance(ids, str):
            ids = [ids]
        
        records = []
        
        try:
            cache_folder = os.path.join(self.cache_dir, database)
            os.makedirs(cache_folder, exist_ok=True)
            for record_id in ids:
                safe_id = record_id.replace("/", "_")
                cache_file = os.path.join(cache_folder, f"{safe_id}.{rettype}.txt")
                if os.path.exists(cache_file):
                    with open(cache_file, "r") as f:
                        text = f.read()
                else:
                    handle = Entrez.efetch(
                        db=database,
                        id=record_id,
                        rettype=rettype,
                        retmode="text"
                    )
                    text = handle.read()
                    handle.close()
                    with open(cache_file, "w") as f:
                        f.write(text)

                if parse and rettype in ["fasta", "gb"]:
                    seq_records = list(SeqIO.parse(io.StringIO(text), rettype))
                    records.extend(seq_records)
                else:
                    records.append(text)

            return records
            
        except Exception as e:
            logger.error(f"Fetch failed: {e}")
            raise
    
    def search_and_fetch(self,
                        query: str,
                        database: str = "nucleotide",
                        retmax: int = 5,
                        field: Optional[str] = None) -> List[SeqRecord]:
        """
        Convenience method to search and fetch in one operation.
        
        Args:
            query: Search query
            database: NCBI database
            retmax: Maximum results
            field: Search field
            
        Returns:
            List of SeqRecord objects
        """
        ids = self.search_database(query, database, retmax, field)
        
        if not ids:
            raise InvalidQueryException(query, database)
        
        return self.fetch_records(ids, database, parse=True)


# Example usage
if __name__ == "__main__":
    # Initialize tools
    seq_tools = SequenceTools()
    
    # Example DNA sequence
    dna_seq = "ATGGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTGA"
    
    # Basic operations
    print("DNA Sequence:", dna_seq)
    print("RNA Sequence:", seq_tools.dna_to_rna(dna_seq))
    print("Reverse Complement:", seq_tools.reverse_complement(dna_seq))
    print(f"GC Content: {seq_tools.gc_content(dna_seq):.2f}%")
    print(f"Melting Temperature: {seq_tools.melting_temperature(dna_seq):.2f}°C")
    print("Protein:", seq_tools.translate(dna_seq))
    
    # Sequence analysis
    analysis = seq_tools.analyze_sequence(dna_seq, SequenceType.DNA)
    print(f"\nSequence Analysis:")
    print(f"  Length: {analysis.length} bp")
    print(f"  GC Content: {analysis.gc_content:.2f}%")
    print(f"  Molecular Weight: {analysis.molecular_weight:.0f} Da")
    
    # Alignment example
    seq1 = "ATCGATCGATCG"
    seq2 = "ATCGTTCGATCG"
    
    result = SequenceAlignment.align_sequences(
        [seq1, seq2], 
        ["seq1", "seq2"],
        seq_type="nucleotide"
    )
    
    print(f"\nAlignment Score: {result.score}")
    print(f"Identity: {result.identity:.2f}%")
    print(f"Gaps: {result.gaps}")
    
    # NCBI example (requires valid email)
    # ncbi = NCBITools("your.email@example.com")
    # records = ncbi.search_and_fetch("BRCA1 human", database="nucleotide", retmax=3)
    # for record in records:
    #     print(f"ID: {record.id}, Description: {record.description[:50]}...")