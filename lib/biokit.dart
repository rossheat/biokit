/// BioKit is a Dart package for object-orientated Bioinformatics.
library biokit;

import 'dart:collection';
import 'package:pdf/pdf.dart';
import 'package:pdf/widgets.dart' as pw;
import 'dart:convert';
import 'dart:io';
import 'dart:math';

/// @nodoc
const String kDNA = 'dna';

/// @nodoc
const String kRNA = 'rna';

/// @nodoc
const String kPep = 'pep';

/// @nodoc
const String kAASeq = 'aaSeq';

/// @nodoc
const String kStartIndex = 'startIndex';

/// @nodoc
const String kEndIndex = 'endIndex';

/// @nodoc
const String kMatchCount = 'matchCount';

/// @nodoc
const String kMatchIndices = 'matchIndices';

/// @nodoc
const String kMatch = 'match';

/// @nodoc
const String kSeq = 'seq';

/// @nodoc
const String kId = 'id';

/// @nodoc
const String kDesc = 'desc';

/// The four DNA nucleotides.
const List<String> dnaNucs = ["A", "T", "G", "C"];

/// The four RNA nucleotides.
const List<String> rnaNucs = ["A", "U", "G", "C"];

/// The three translation stop codons.
const List<String> translationStopCodons = ["UGA", "UAA", "UAG"];

/// The twenty amino acids.
const List<String> aminoAcids = [
  "F",
  "S",
  "Y",
  "C",
  "W",
  "L",
  "P",
  "H",
  "Q",
  "I",
  "M",
  "T",
  "N",
  "K",
  "R",
  "V",
  "A",
  "D",
  "E",
  "G",
  "X" // Stop codon
];

/// DNA transitions.
const List<String> dnaTransitions = ['CT', 'TC', 'AG', 'GA'];

/// DNA transversions.
const List<String> dnaTransversions = ['GT', 'TG', 'AC', 'CA', 'AT', 'TA', 'GC', 'CG'];

/// Complementary DNA nucleotides.
const Map<String, String> dnaCompNucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'};

/// Complementary RNA nucelotides.
const Map<String, String> rnaCompNucs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'};

/// DNA codons to amino acids.
///
/// `"X"` represents the occurance of a stop codon.
const Map<String, String> dnaCodonToAA = {
  "TTT": "F",
  "TTC": "F",
  "TTA": "L",
  "TTG": "L",
  "TCT": "S",
  "TCC": "S",
  "TCA": "S",
  "TCG": "S",
  "TAT": "Y",
  "TAC": "Y",
  "TAA": "X",
  "TAG": "X",
  "TGT": "C",
  "TGC": "C",
  "TGA": "X",
  "TGG": "W",
  "CTT": "L",
  "CTC": "L",
  "CTA": "L",
  "CTG": "L",
  "CCT": "P",
  "CCC": "P",
  "CCA": "P",
  "CCG": "P",
  "CAT": "H",
  "CAC": "H",
  "CAA": "Q",
  "CAG": "Q",
  "CGT": "R",
  "CGC": "R",
  "CGA": "R",
  "CGG": "R",
  "ATT": "I",
  "ATC": "I",
  "ATA": "I",
  "ATG": "M",
  "ACT": "T",
  "ACC": "T",
  "ACA": "T",
  "ACG": "T",
  "AAT": "N",
  "AAC": "N",
  "AAA": "K",
  "AAG": "K",
  "AGT": "S",
  "AGC": "S",
  "AGA": "R",
  "AGG": "R",
  "GTT": "V",
  "GTC": "V",
  "GTA": "V",
  "GTG": "V",
  "GCT": "A",
  "GCC": "A",
  "GCA": "A",
  "GCG": "A",
  "GAT": "D",
  "GAC": "D",
  "GAA": "E",
  "GAG": "E",
  "GGT": "G",
  "GGC": "G",
  "GGA": "G",
  "GGG": "G"
};

/// RNA codons to amino acids.
///
/// `"X"` represents the occurance of a stop codon.
const Map<String, String> rnaCodonToAA = {
  "UUU": "F",
  "UUC": "F",
  "UUA": "L",
  "UUG": "L",
  "UCU": "S",
  "UCC": "S",
  "UCA": "S",
  "UCG": "S",
  "UAU": "Y",
  "UAC": "Y",
  "UAA": "X",
  "UAG": "X",
  "UGU": "C",
  "UGC": "C",
  "UGA": "X",
  "UGG": "W",
  "CUU": "L",
  "CUC": "L",
  "CUA": "L",
  "CUG": "L",
  "CCU": "P",
  "CCC": "P",
  "CCA": "P",
  "CCG": "P",
  "CAU": "H",
  "CAC": "H",
  "CAA": "Q",
  "CAG": "Q",
  "CGU": "R",
  "CGC": "R",
  "CGA": "R",
  "CGG": "R",
  "AUU": "I",
  "AUC": "I",
  "AUA": "I",
  "AUG": "M",
  "ACU": "T",
  "ACC": "T",
  "ACA": "T",
  "ACG": "T",
  "AAU": "N",
  "AAC": "N",
  "AAA": "K",
  "AAG": "K",
  "AGU": "S",
  "AGC": "S",
  "AGA": "R",
  "AGG": "R",
  "GUU": "V",
  "GUC": "V",
  "GUA": "V",
  "GUG": "V",
  "GCU": "A",
  "GCC": "A",
  "GCA": "A",
  "GCG": "A",
  "GAU": "D",
  "GAC": "D",
  "GAA": "E",
  "GAG": "E",
  "GGU": "G",
  "GGC": "G",
  "GGA": "G",
  "GGG": "G"
};

/// Amino acids and their monoisotopic mass.
const Map<String, double> aaToMonoMass = {
  "A": 71.03711,
  "C": 103.00919,
  "D": 115.02694,
  "E": 129.04259,
  "F": 147.06841,
  "G": 57.02146,
  "H": 137.05891,
  "I": 113.08406,
  "K": 128.09496,
  "L": 113.08406,
  "M": 131.04049,
  "N": 114.04293,
  "P": 97.05276,
  "Q": 128.05858,
  "R": 156.10111,
  "S": 87.03203,
  "T": 101.04768,
  "V": 99.06841,
  "W": 186.07931,
  "Y": 163.06333,
  "X": 0,
};

/// A collection of error functions.
///
/// @nodoc
class Errors {
  /// Throws a sequence error.
  static void invalidSeq({required String mon, required int idx, required String type}) {
    throw ("Invalid ${type.toUpperCase()} Sequence Error. Character '$mon' found at index position $idx (zero-based) is not a valid ${type.toUpperCase()} monomer.");
  }

  /// Throws a sequence type error.
  static void invalidType({required String type}) {
    throw ("Invalid Sequence Type Error. '$type' is not a valid sequence type. Please enter the argument 'dna', 'rna' or 'pep' for [type].");
  }

  /// Throws an invalid FASTA arguement error..
  static void invalidFASTAArgs() {
    throw ('Invalid readFASTA arguement specification. Please specify either a [path] or a [str] in which to read FASTA formatted data from.');
  }

  /// Throws an unequal sequence type error.
  static void unequalSeqTypes({required String func}) {
    throw ('$func function requires sequences must be of the same type.');
  }

  /// Throws an unequal sequence length error.
  static void unequalSeqLens({required String func}) {
    throw ('$func function requires sequences must be of the same length.');
  }

  /// Throws an invalid sequence length error.
  static void invalidSeqLen({required String param, required int min, int? max}) {
    if (max != null) {
      throw ('$param must have between $min and $max characters.');
    } else {
      throw ('$param must have $min or more characters.');
    }
  }
}

/// A collection of helper functions for bioinformatics.
class Utils {
  /// Returns a protein sequence in FASTA format using it's [uniprotId].
  ///
  /// Requires network access.
  static Future<String> uniprotIdToFASTA({required String uniprotId}) async {
    Uri uri = Uri.parse('http://www.uniprot.org/uniprot/$uniprotId.fasta');

    var request = await HttpClient().getUrl(uri);
    var response = await request.close();

    await for (var contents in response.transform(Utf8Decoder())) {
      return contents;
    }
    return 'Error retrieving protein with uniprot ID $uniprotId';
  }

  /// Returns sequences from either a file [path] or [str] in FASTA format.
  ///
  /// Each map contains a sequence `seq`, ID `id`, and description `desc`.
  static Future<List<Map<String, String>>> readFASTA({String? path, String? str}) async {
    String? contents;
    if (path != null) {
      contents = await File(path).readAsString();
    } else if (str != null) {
      contents = str;
    } else {
      Errors.invalidFASTAArgs();
    }

    List<String> lines = contents!.split('\n');
    int seqCount = 0;

    List<Map<String, String>> fastaMaps = [];
    Map<String, String> currentMap = {};

    for (var line in lines) {
      if (line.startsWith('>')) {
        if (seqCount != 0) {
          fastaMaps.add(currentMap);
          currentMap = {};
        }
        seqCount++;
        currentMap[kSeq] = '';

        String topLine = line.split('>')[1];
        List<String> topLineList = topLine.split(' ');

        currentMap[kId] = topLineList.first;
        currentMap[kDesc] = topLineList.sublist(1, topLineList.length).join();
      } else {
        currentMap[kSeq] = currentMap[kSeq]! + line;
      }
    }
    fastaMaps.add(currentMap);
    return fastaMaps;
  }

  /// Returns a regex valid version of a biological [motif].
  ///
  /// ```dart
  /// motifToRe(motif: 'N{P}[ST]{P}') == 'N[^P][S|T|][^P]'
  /// ```
  static String motifToRe({required String motif}) {
    String re = '';

    List<String> chars = motif.split('');

    bool inBrac = false;
    for (String char in chars) {
      if (char == '[') {
        inBrac = true;
        re += char;
      } else if (char == ']') {
        inBrac = false;
        re += char;
      } else {
        if (inBrac) {
          re += char + '|';
        } else {
          re += char;
        }
      }
    }
    return re.replaceAll('{', '[^').replaceAll('}', ']');
  }
}

/// A model representation of a biological sequence.
class Sequence {
  /// The sequence.
  late final String _seq;

  /// The sequence type.
  late final String _type;

  /// The number of monomers.
  late final int _len;

  /// The name.
  late String _name;

  /// The ID.
  late String _id;

  /// The description.
  late String _desc;

  /// Creates a `Sequence` object.
  Sequence._({required String seq, required String type}) {
    this._type = _validateType(type: type);
    this._seq = _validateSeq(seq: seq);
    this._name = 'Default name';
    this._id = 'Default ID';
    this._desc = 'Default description';
  }

  /// Validates and returns a sequence [type].
  String _validateType({required String type}) {
    String lType = type.toLowerCase();
    // Checks if the lowercased type is present in the list
    if (![kDNA, kRNA, kPep].contains(lType)) {
      Errors.invalidType(type: lType);
    }
    return lType;
  }

  /// Validates and returns a [seq].
  String _validateSeq({required String seq}) {
    int seqLen = seq.length;

    if (this is Peptide) {
      // Peptides must be at least one monomer long.
      if (seqLen < 1) {
        Errors.invalidSeqLen(param: 'Peptide', min: 1);
      }
    } else {
      // Nucleotides must at least 6 monomers long.
      if (seqLen < 6) {
        Errors.invalidSeqLen(param: 'Nucelotides', min: 6);
      }
    }

    String uSeq = seq.toUpperCase();

    // Enumerate over each element and confirm it is valid for its type.
    uSeq.split('').asMap().forEach((idx, mon) {
      if (this._type == kDNA) {
        if (!dnaNucs.contains(mon)) {
          Errors.invalidSeq(mon: mon, idx: idx, type: this._type);
        }
      } else if (this.type == kRNA) {
        if (!rnaNucs.contains(mon)) {
          Errors.invalidSeq(mon: mon, idx: idx, type: this._type);
        }
      } else {
        if (!aminoAcids.contains(mon)) {
          Errors.invalidSeq(mon: mon, idx: idx, type: this._type);
        }
      }
    });
    this._len = seqLen;
    return uSeq;
  }

  /// Returns the `String` sequence combination between this and [oSeq].
  String operator +(Sequence oSeq) {
    if (this._type != oSeq.type) {
      Errors.unequalSeqTypes(func: '+');
    }
    return this._seq + oSeq.seq;
  }

  /// Returns information about the this sequence.
  Map<String, dynamic> info() {
    return {
      'seq': seq,
      'type': type,
      'monomers': _len,
      'name': this._name,
      'id': this._id,
      'desc': this._desc
    };
  }

  /// Returns a `String` representation of the [info()] function.
  @override
  String toString() => info().toString();

  /// Returns the sequence.
  String get seq => this._seq;

  /// Returns the sequence type.
  String get type => this._type;

  /// Returns the number of monomers in the sequence.
  int get len => this._len;

  /// Returns the sequence name.
  String get name => this._name;

  /// Returns the sequence ID.
  String get id => this._id;

  /// Returns the sequence description.
  String get desc => this._desc;

  /// Sets the sequence name to [newName].
  set name(String newName) {
    if (newName.length < 2 || 25 > newName.length) {
      Errors.invalidSeqLen(param: 'Name', min: 2, max: 25);
    }
    this._name = newName;
  }

  /// Sets the sequence ID to [newId].
  set id(String newId) {
    if (newId.length < 2 || 30 > newId.length) {
      Errors.invalidSeqLen(param: 'ID', min: 2, max: 30);
    }
    this._id = newId;
  }

  /// Sets the sequence description to [newDesc].
  set desc(String newDesc) {
    if (newDesc.length < 5 || 100 > newDesc.length) {
      Errors.invalidSeqLen(param: 'Description', min: 5, max: 100);
    }
    this._desc = newDesc;
  }

  /// Returns the frequency of each monomer.
  ///
  /// Return the percentage total of each monomer by setting [norm] to `true`.
  /// Inlcude the X 'amino acid' in the count by setting [ignoreStopCodon] to `false`.
  Map<String, double> freq({bool norm = false, bool ignoreStopCodon = true}) {
    Map<String, double> freqMap = {};
    this._seq.split('').forEach((mon) {
      // If this is a peptide and [ignoreStopCodon] is `true`, we want
      // to perform the calulation without taking the stop codon product
      // X into consideration.
      if (this is Peptide && ignoreStopCodon) {
        if (mon != 'X') {
          if (freqMap.containsKey(mon)) {
            freqMap[mon] = freqMap[mon]! + 1;
          } else {
            freqMap[mon] = 1;
          }
        }
        // If either this is not a Peptide or [ignoreStopCodon] is `false`,
        // then proceed without checking if the character is X.
      } else {
        if (freqMap.containsKey(mon)) {
          freqMap[mon] = freqMap[mon]! + 1;
        } else {
          freqMap[mon] = 1;
        }
      }
    });

    if (this is Peptide && ignoreStopCodon) {
      if (norm) {
        // The sequence is a Peptide and [ignoreStopCodon] is `true` and [norm] is `true`.
        // Return normalised results based on length without the X character.
        return freqMap.map(
          (key, value) => MapEntry(
            key,
            // Use the length of the peptide in absence of the X character.
            double.parse(((value / lenMinus(monomer: 'X')) * 100).toStringAsFixed(1)),
          ),
        );
      }
      // The sequence is not a Peptide or [ignoreStopCodon] is `false`.
    } else {
      if (norm) {
        // Return normalised results based on sequence length.
        return freqMap.map(
          (key, value) => MapEntry(
            key,
            double.parse(((value / this._len) * 100).toStringAsFixed(1)),
          ),
        );
      }
    }
    // The sequence is not a Peptide or [ignoreStopCodon] is `false`, or [norm] is `false`.
    // Return non-normalised results based on sequence length.
    return freqMap;
  }

  /// Returns the reversed version of this sequence.
  String _reversed({required String seq}) => seq.split('').reversed.join('');

  /// Returns the reversed version of this sequence.
  String reverse() => this._seq.split('').reversed.join('');

  /// Returns every possible combination of this sequence.
  ///
  /// Sort the result from longest to shortest by setting [sorted] to `true`.
  ///
  /// ```dart
  /// Peptide pep = Peptide(seq: 'ATG');
  /// pep.combinations()
  /// // > [A, AT, ATG, T, TG];
  /// ```
  List<String> combinations({sorted = false}) {
    List<String> listSeq = this._seq.split("");
    List<String> combinations = [];

    // Generates all possible combinations of this sequence.
    for (int i = 0; i < listSeq.length; i++) {
      if (i != listSeq.length - 1) {
        combinations.add(listSeq[i]);
      }
      List<String> temp = [listSeq[i]];
      for (int j = i + 1; j < listSeq.length; j++) {
        temp.add(listSeq[j]);
        combinations.add(temp.join());
      }
    }
    if (sorted) {
      // Sorts with longest combination first.
      combinations.sort((b, a) => a.length.compareTo(b.length));
      return combinations;
    }
    return combinations;
  }

  /// Returns the indices of all [motif] matches.
  ///
  /// Prevent overlapping matches by setting [overlap] to `false`.
  Map<String, dynamic> findMotif({required String motif, overlap = true}) {
    List<Map<String, dynamic>> matchData = [];
    Map<String, dynamic> matchMotifMap = {};

    // 'ATGC' -> 'ATGC', 'N{P}[ST]{P}' -> 'N[^P][S|T|][^P]'
    String tempRegexMotif = Utils.motifToRe(motif: motif);

    // ?= is the regex lookahead operator and allows for matches to overlap.
    RegExp regexMotif = overlap == true ? RegExp('(?=$tempRegexMotif)') : RegExp(tempRegexMotif);

    Iterable<RegExpMatch> allMatches = regexMotif.allMatches(seq);
    for (RegExpMatch match in allMatches) {
      matchData.add({
        kMatch: motif,
        kStartIndex: match.start,
        kEndIndex: match.start + motif.length - 1,
      });
    }
    matchMotifMap[kMatchCount] = allMatches.length;
    matchMotifMap[kMatchIndices] = matchData;
    return matchMotifMap;
  }

  /// Returns the number of monomer differences between this sequence and [oSeq].
  int difference({required Sequence oSeq}) {
    if (this._len != oSeq.len) {
      Errors.unequalSeqLens(func: 'difference');
    }
    if (this._type != oSeq.type) {
      Errors.unequalSeqTypes(func: 'difference');
    }

    int differenceCount = 0;
    this._seq.split('').asMap().forEach((idx, mon) {
      if (mon != oSeq.seq[idx]) {
        differenceCount++;
      }
    });
    return differenceCount;
  }

  /// Returns this sequence with all occurrences of [motif] removed.
  String splice({required String motif}) {
    String vMotif = _validateSeq(seq: motif);
    return seq.replaceAll(vMotif, '');
  }

  /// Returns the longest shared motif between this sequence and [oSeq].
  String sharedMotif({required Sequence oSeq}) {
    if (this._type != oSeq.type) {
      Errors.unequalSeqTypes(func: 'sharedMotif');
    }
    // Return all possible combinates of this sequence sorted from longest to shortest.
    List<String> combos = combinations(sorted: true);

    String longestShared = '';

    for (var comb in combos) {
      // Assume the combination is contained in the other sequence unless set otherwise.
      bool isMatch = true;
      if (!oSeq.seq.contains(comb)) {
        // The combination is not contained in the other sequence.
        isMatch = false;
        break;
      }
      if (isMatch) {
        // The combination must be present in the other sequence.
        // Because the combinations were sorted from longest to
        // shortest, this must be the longest shared motif.
        longestShared = comb;
        break;
      }
    }
    return longestShared;
  }

  /// Returns a number normalized by a [total] to a number of decimal [places].
  double _norm({required int number, required int total, required int places}) =>
      double.parse(((number / total) * 100).toStringAsFixed(1));

  /// Returns the length of this sequence with [monomer] removed.
  int lenMinus({required String monomer}) {
    int newLen = 0;
    for (var mon in this._seq.split('')) {
      if (mon != monomer) {
        newLen++;
      }
    }
    return newLen;
  }

  /// Saves this sequence to a local file in FASTA format.
  ///
  /// ```dart
  /// dna.toFASTA(outputPath: '../deliverables', filename: 'fasta_output');
  /// ```
  Future<void> toFASTA({required String outputPath, required String filename}) async {
    final String metadata = '>${this._id} ${this._desc}\n';
    String seq = '';
    int lineLen = 0;
    for (var mon in this._seq.split('')) {
      lineLen++;
      // Each sequence line should be at most 60 characters long.
      if (lineLen % 60 == 0) {
        // [lineLen] is a multiple of 60 so start a new line.
        seq += '$mon\n';
      } else {
        seq += '$mon';
      }
    }
    final String fileData = metadata + seq;
    final String filen = filename.toLowerCase().trim();
    await File('$outputPath/$filen.txt').writeAsString(fileData);
  }
}

/// A model representation of a nucleotide sequence.
class Nucleotides extends Sequence {
  /// Creates a `Nucleotides` object.
  Nucleotides._({required String seq, required String type}) : super._(seq: seq, type: type);

  /// Returns the translated version of this sequence.
  ///
  /// Return the reverse complementary strand by setting [revComp] to `true`.
  /// Alter the starting indexing by setting [startIdx].
  Map<String, dynamic> translate({revComp = false, startIdx = 0}) {
    // If [revComp] is `true`, return the reverse complementary sequence to this sequence.
    String seq = revComp ? complementary(rev: true) : this.seq;

    String aaSeq = '';
    // Loop through the sequence in batches of three characters at a time.
    for (var i = startIdx; i < seq.length - 2; i += 3) {
      // Grab three characters from the sequence.
      String codon = seq.substring(i, i + 3);
      // Retrive the matching amino acid from the correct codon table based on this sequence type.
      aaSeq += this.type == kDNA ? dnaCodonToAA[codon]! : rnaCodonToAA[codon]!;
    }
    return {kAASeq: aaSeq, 'nucCount': seq.length - startIdx - 1, 'aaCount': aaSeq.length};
  }

  /// Returns the frequency of a specified [codon].
  ///
  /// Scans in batches of three monomers per step.
  /// The exact [codon] must be present in a batch to be detected.
  int codonFreq({required String codon}) {
    int codonFreq = 0;
    // Loop through the sequence in batches of three characters at a time.
    for (var i = 0; i < this._len - 2; i += 3) {
      // Grab three characters from the sequence and check if they match [codon].
      String fetchedCodon = this._seq.substring(i, i + 3);
      if (codon == fetchedCodon) {
        codonFreq++;
      }
    }
    return codonFreq;
  }

  /// Returns the complementary strand to this sequence.
  ///
  /// Return the reversed complementary strand by setting [rev] to `true`.
  String complementary({bool rev = false}) {
    String compSeq = this
        .seq
        .split('')
        .map((nuc) => this.type == kDNA ? dnaCompNucs[nuc] : rnaCompNucs[nuc])
        .join();
    return rev ? super._reversed(seq: compSeq) : compSeq;
  }

  /// Returns the percentage of Guanine and Cytosine nucleotides in this sequence.
  double gcContent() {
    int gcCount = 0;
    this.seq.split('').forEach((nuc) {
      if (nuc == 'G' || nuc == 'C') {
        gcCount++;
      }
    });
    return double.parse(
      ((gcCount / this._len) * 100).toStringAsFixed(2),
    );
  }

  /// Returns six reading frames from this sequence.
  List<String> readingFrames() {
    List<String> readingFrames = [];

    for (var i = 0; i < 3; i++) {
      // Generate three forwards reading frames.
      readingFrames.add(translate(startIdx: i)[kAASeq]);
      // Generate three backwards reading frames.
      readingFrames.add(translate(revComp: true, startIdx: i)[kAASeq]);
    }
    return readingFrames;
  }

  /// Returns protein sequences from a single [aaSeq] sequence.
  ///
  /// Commonly used after generating [readingFrames()].
  List<String> readingFrameToProteins({required String aaSeq}) {
    List<String> currentProtein = [];
    List<String> proteins = [];

    for (var aa in aaSeq.split('')) {
      if (aa == 'X') {
        if (currentProtein != []) {
          for (var pro in currentProtein) {
            proteins.add(pro);
          }
          currentProtein = [];
        }
      } else {
        if (aa == 'M') {
          currentProtein.add('');
        }
        for (var i = 0; i < currentProtein.length; i++) {
          currentProtein[i] += aa;
        }
      }
    }
    return proteins;
  }

  /// Returns proteins from this sequence.
  List<String> proteins({bool unique = false}) {
    // Generate siz reading frames from this sequence.
    List<String> frames = readingFrames();
    List<String> proteins = [];
    // Loop through each of the six reading frames.
    for (var frame in frames) {
      // Translate the reading frame into protein(s).
      List<String> tempProteins = readingFrameToProteins(aaSeq: frame);
      if (proteins != []) {
        // If atleast one protein has been generated, at it/them to [proteins].
        proteins.addAll(tempProteins);
      }
    }
    // sort the proteins from longest to shortest.
    proteins.sort((b, a) => a.length.compareTo(b.length));
    if (unique) {
      // If [unique] is `true` only return the unique proteins.
      return proteins.toSet().toList();
    }
    return proteins;
  }

  /// Returns the molecular weight (kDA) of this sequence.
  double molWeight() => this._len * 0.33;
}

/// A model representation of a DNA sequence.
class DNA extends Nucleotides {
  /// Creates a `DNA` object.
  DNA({required String seq}) : super._(seq: seq, type: 'dna');

  /// Returns the molecular weight (kDA) of this sequence if it were in double helical form.
  double doubleHelixMolWeight() => molWeight() * 2;

  /// Returns the number of turns in this sequence if it were in double helical form.
  double doubleHexlixTurns() => this._len / 10;

  /// Returns the length (nm) of this sequence if it were in double helical form.
  double doubleHelixGeoLen() => this._len * 0.34;

  /// Returns the transcribed version of this sequence.
  String transcribe() => this.seq.replaceAll('T', 'U');

  /// Returns the restriction sites of this sequence.
  ///
  /// Alter the length of the restriction sites by modifying [minSiteLen] and [maxSiteLen].
  Map<String, List<Map<String, int>>> restrictionSites({
    int minSiteLen = 4,
    int maxSiteLen = 8,
  }) {
    // Fetch a reverse complementary sequence to this sequence.
    String revComp = complementary(rev: true);

    // Fetch all possible combinations of this sequence.
    List<String> origSeqCombos = combinations();

    // Fetch only the combinations whose length is equato to or between [minSiteLen] and [maxSiteLen].
    Iterable<String> seqIter =
        origSeqCombos.where((seq) => (seq.length >= minSiteLen) && (seq.length <= maxSiteLen));

    // Fecth only the unique combinations.
    List<String> restSeqs = seqIter.toSet().toList();

    Map<String, List<Map<String, int>>> restSiteSeqs = {};
    List<Map<String, int>> restSiteLocations = [];

    // Loop through each combination.
    for (String restSeq in restSeqs) {
      restSiteLocations = [];
      RegExp regRestExp = RegExp(restSeq);
      // Find all matches between the combination and the reverse complementary sequence.
      Iterable<RegExpMatch> matches = regRestExp.allMatches(revComp);
      // Loop through each match.
      for (RegExpMatch match in matches) {
        int startIdx = this.seq.length - (restSeq.length + match.start);
        int endIdx = startIdx + restSeq.length;
        // Detect palindrome.
        // Check if the sub-sequence matches the current combination.
        if (this.seq.substring(startIdx, endIdx) == restSeq) {
          restSiteLocations.add({kStartIndex: startIdx, kEndIndex: endIdx});
          restSiteSeqs[restSeq] = restSiteLocations;
        }
        ;
      }
    }
    return restSiteSeqs;
  }

  /// Returns the transition-transversion ratio between this sequence and [oSeq].
  double tranRatio({required DNA oSeq}) {
    if (this._len != oSeq._len) {
      Errors.unequalSeqLens(func: 'tranRatio');
    }

    int transitionCount = 0;
    int transversionCount = 0;

    this.seq.split('').asMap().forEach((idx, nuc) {
      if (nuc != oSeq.seq[idx]) {
        // 'A' + 'T' => 'AT'
        if (dnaTransitions.contains(nuc + oSeq.seq[idx])) {
          transitionCount++;
        } else if (dnaTransversions.contains(nuc + oSeq.seq[idx])) {
          transversionCount++;
        }
      }
    });
    // Cant divide by 0 (inf).
    return transversionCount == 0 ? -1 : (transitionCount / transversionCount);
  }

  /// Returns a `DNA` object with a specified length of [len].
  static DNA random({required int len}) {
    Random _rand = Random();
    String dnaNucsStr = dnaNucs.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => dnaNucsStr.codeUnitAt(
          _rand.nextInt(dnaNucsStr.length),
        ),
      ),
    );
    return DNA(seq: seq);
  }

  /// Generates and saves a DNA analysis report to [outputPath].
  ///
  /// Add your name to the report by setting [creatorName].
  /// Add a title to the report  by setting [reportTitle].
  ///
  /// ```dart
  /// DNA dna = DNA(seq: 'ATGCGA');
  /// dna.report(outputPath: '../deliverables', reportTitle: 'My Report', creatorName: 'John Doe');
  /// ```
  Future<void> report(
      {required String outputPath,
      required String creatorName,
      required String reportTitle}) async {
    _genReport(outputPath, creatorName, reportTitle);
  }

  Future<void> _genReport(outputPath, creatorName, reportTitle) async {
    final DateTime now = new DateTime.now();
    final DateTime date = new DateTime(now.year, now.month, now.day);

    final int monomers = this._len;
    final double gcCon = gcContent();
    final double mWeight = molWeight();
    final double turns = doubleHexlixTurns();
    final double geoLength = doubleHelixGeoLen();
    final Map nucFreqCount = freq();
    final Map nucFreqPerc = freq(norm: true);

    // DNA
    const dnaRightTableHeaders = ['Nucleotide', 'Count', 'Percent Total (%)'];
    var dnaRightTableData = [
      ['A', nucFreqCount['A'].toStringAsFixed(0), nucFreqPerc['A']],
      ['T', nucFreqCount['T'].toStringAsFixed(0), nucFreqPerc['T']],
      ['G', nucFreqCount['G'].toStringAsFixed(0), nucFreqPerc['G']],
      ['C', nucFreqCount['C'].toStringAsFixed(0), nucFreqPerc['C']],
    ];
    var dnaTopTableHeaders = [
      'Nucelotides',
      'GC Content (%)',
      'mWeight (kDa)',
      'DHelix Length (nm)',
      'DHelix Turns'
    ];
    var dnaTopTableData = [
      monomers,
      gcCon,
      mWeight.toStringAsFixed(1),
      geoLength.toStringAsFixed(1),
      turns.toStringAsFixed(1),
    ];

    // RNA
    RNA rna = RNA(seq: transcribe());

    final int numOfCodons = (this._len / 3).floor();
    final int aug = rna.codonFreq(codon: 'AUG');
    final int uag = rna.codonFreq(codon: 'UAG');
    final int uaa = rna.codonFreq(codon: 'UAA');
    final int uga = rna.codonFreq(codon: 'UGA');
    final int other = numOfCodons - aug - uag - uaa - uga;

    var rnaTopTableHeaders = [
      'Codons',
      'Total Start Codons',
      'Total Stop Codons',
    ];
    var rnaTopTableData = [
      numOfCodons,
      aug,
      uag + uaa + uga,
    ];

    const rnaRightTableHeaders = ['Codon', 'Count', 'Percent Total (%)'];
    var rnaRightTableData = [
      ['Other', other, _norm(number: other, total: numOfCodons, places: 1)],
      ['AUG', aug, _norm(number: aug, total: numOfCodons, places: 1)],
      ['UAG', uag, _norm(number: uag, total: numOfCodons, places: 1)],
      ['UAA', uaa, _norm(number: uaa, total: numOfCodons, places: 1)],
      ['UGA', uga, _norm(number: uga, total: numOfCodons, places: 1)],
    ];

    // Peptide
    Peptide pep = Peptide(seq: translate()[kAASeq]);

    var peptideTopTableHeaders = [
      'Amino Acids',
      'Monoisotopic Mass (kDa)',
    ];
    var peptideTopTableData = [
      pep.lenMinus(monomer: 'X'),
      pep.monoMass(roundTo: 1, kDa: true),
    ];

    Map<String, double> aaFreq = pep.freq();
    Map<String, double> aaFreqPct = pep.freq(norm: true);

    var sortedKeys = aaFreq.keys.toList(growable: false)
      ..sort((k2, k1) => aaFreq[k1]!.compareTo(aaFreq[k2]!));
    LinkedHashMap sortedMap =
        new LinkedHashMap.fromIterable(sortedKeys, key: (k) => k, value: (k) => aaFreq[k]);

    const pepRightTableHeaders = ['Top 5 AAs', 'Count', 'Percent Total (%)'];
    List<List<dynamic>> pepRightTableData = [];
    for (var key in sortedMap.keys.take(5)) {
      pepRightTableData.add([
        key,
        aaFreq[key]!.toInt(),
        aaFreqPct[key],
      ]);
    }

    final baseColor = PdfColors.cyan;
    final document = pw.Document();

    document.addPage(
      pw.Page(
        pageFormat: null,
        build: (context) {
          final dnaBarChart = pw.Chart(
            left: pw.Container(
              alignment: pw.Alignment.topCenter,
              margin: const pw.EdgeInsets.only(right: 5, top: 10),
              child: pw.Transform.rotateBox(
                angle: pi / 2,
                child: pw.Text('Percent Total (%)'),
              ),
            ),
            grid: pw.CartesianGrid(
              xAxis: pw.FixedAxis.fromStrings(
                List<String>.generate(
                    dnaRightTableData.length,
                    // ignore: avoid_as
                    (index) => dnaRightTableData[index][0] as String),
                marginStart: 30,
                marginEnd: 30,
                ticks: true,
              ),
              yAxis: pw.FixedAxis(
                [0, 20, 40, 60, 80, 100],
                format: (v) => '$v\%',
                divisions: true,
              ),
            ),
            datasets: [
              pw.BarDataSet(
                color: baseColor,
                legend: dnaRightTableHeaders[2],
                width: 30,
                offset: 0,
                borderColor: baseColor,
                data: List<pw.LineChartValue>.generate(
                  dnaRightTableData.length,
                  (i) {
                    // ignore: avoid_as
                    final v = dnaRightTableData[i][2] as num;
                    return pw.LineChartValue(i.toDouble(), v.toDouble());
                  },
                ),
              ),
            ],
          );

          final dnaRightTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            border: null,
            headers: dnaRightTableHeaders,
            data: dnaRightTableData,
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          final dnaTopTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
              3: pw.Alignment.center,
              4: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
              3: pw.Alignment.center,
              4: pw.Alignment.center,
            },
            border: null,
            headers: dnaTopTableHeaders,
            data: [dnaTopTableData],
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          final rnaTopTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
              3: pw.Alignment.center,
              4: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
              3: pw.Alignment.center,
              4: pw.Alignment.center,
            },
            border: null,
            headers: rnaTopTableHeaders,
            data: [rnaTopTableData],
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          final rnaBarChart = pw.Chart(
            left: pw.Container(
              alignment: pw.Alignment.topCenter,
              margin: const pw.EdgeInsets.only(right: 5, top: 10),
              child: pw.Transform.rotateBox(
                angle: pi / 2,
                child: pw.Text('Percent Total (%)'),
              ),
            ),
            grid: pw.CartesianGrid(
              xAxis: pw.FixedAxis.fromStrings(
                List<String>.generate(
                    rnaRightTableData.length,
                    // ignore: avoid_as
                    (index) => rnaRightTableData[index][0] as String),
                marginStart: 30,
                marginEnd: 30,
                ticks: true,
              ),
              yAxis: pw.FixedAxis(
                [0, 20, 40, 60, 80, 100],
                format: (v) => '$v\%',
                divisions: true,
              ),
            ),
            datasets: [
              pw.BarDataSet(
                color: baseColor,
                legend: rnaRightTableHeaders[2],
                width: 30,
                offset: 0,
                borderColor: baseColor,
                data: List<pw.LineChartValue>.generate(
                  rnaRightTableData.length,
                  (i) {
                    // ignore: avoid_as
                    final v = rnaRightTableData[i][2] as num;
                    return pw.LineChartValue(i.toDouble(), v.toDouble());
                  },
                ),
              ),
            ],
          );

          final rnaRightTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            border: null,
            headers: rnaRightTableHeaders,
            data: rnaRightTableData,
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          final pepTopTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
            },
            border: null,
            headers: peptideTopTableHeaders,
            data: [peptideTopTableData],
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          final pepBarChart = pw.Chart(
            left: pw.Container(
              alignment: pw.Alignment.topCenter,
              margin: const pw.EdgeInsets.only(right: 5, top: 10),
              child: pw.Transform.rotateBox(
                angle: pi / 2,
                child: pw.Text('Percent Total (%)'),
              ),
            ),
            grid: pw.CartesianGrid(
              xAxis: pw.FixedAxis.fromStrings(
                List<String>.generate(
                    pepRightTableData.length,
                    // ignore: avoid_as
                    (index) => pepRightTableData[index][0] as String),
                marginStart: 30,
                marginEnd: 30,
                ticks: true,
              ),
              yAxis: pw.FixedAxis(
                [0, 20, 40, 60, 80, 100],
                format: (v) => '$v\%',
                divisions: true,
              ),
            ),
            datasets: [
              pw.BarDataSet(
                color: baseColor,
                legend: pepRightTableHeaders[2],
                width: 25,
                offset: 0,
                borderColor: baseColor,
                data: List<pw.LineChartValue>.generate(
                  pepRightTableData.length,
                  (i) {
                    // ignore: avoid_as
                    final v = pepRightTableData[i][2] as num;
                    return pw.LineChartValue(i.toDouble(), v.toDouble());
                  },
                ),
              ),
            ],
          );

          final pepRightTable = pw.Table.fromTextArray(
            headerAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            cellAlignments: {
              0: pw.Alignment.center,
              1: pw.Alignment.center,
              2: pw.Alignment.center,
            },
            border: null,
            headers: pepRightTableHeaders,
            data: pepRightTableData,
            headerStyle: pw.TextStyle(
              color: PdfColors.white,
              fontWeight: pw.FontWeight.bold,
            ),
            headerDecoration: pw.BoxDecoration(
              color: baseColor,
            ),
            rowDecoration: pw.BoxDecoration(
              border: pw.Border(
                bottom: pw.BorderSide(
                  color: baseColor,
                  width: .5,
                ),
              ),
            ),
          );

          return pw.Column(
            children: [
              pw.Text(
                reportTitle,
                style: pw.TextStyle(
                  fontSize: 25,
                ),
              ),
              pw.SizedBox(height: 2),
              pw.Text(
                'by $creatorName on ${date.day}/${date.month}/${date.year} with BioKit.org',
                style: pw.TextStyle(color: PdfColors.grey, fontSize: 11),
              ),
              pw.SizedBox(height: 15),
              pw.Row(children: [
                pw.Text(id),
                pw.Text(name),
                pw.Text('${this.seq.substring(0, 3)} ... ${seq.substring(this.len - 3, this.len)}')
              ], mainAxisAlignment: pw.MainAxisAlignment.spaceBetween),
              // DNA
              pw.SizedBox(height: 15),
              dnaTopTable,
              pw.Row(
                crossAxisAlignment: pw.CrossAxisAlignment.start,
                children: [
                  pw.Container(
                      child: pw.Expanded(flex: 1, child: dnaBarChart), height: 136, width: 240),
                  pw.SizedBox(width: 10),
                  pw.Expanded(flex: 1, child: dnaRightTable),
                ],
              ),
              //RNA
              pw.SizedBox(height: 20),
              rnaTopTable,
              pw.Row(
                crossAxisAlignment: pw.CrossAxisAlignment.start,
                children: [
                  pw.Container(
                      child: pw.Expanded(flex: 1, child: rnaBarChart), height: 160, width: 270),
                  pw.SizedBox(width: 10),
                  pw.Expanded(flex: 1, child: rnaRightTable),
                ],
              ),
              //Peptide
              pw.SizedBox(height: 20),
              pepTopTable,
              pw.Row(
                crossAxisAlignment: pw.CrossAxisAlignment.start,
                children: [
                  pw.Container(
                      child: pw.Expanded(flex: 1, child: pepBarChart), height: 160, width: 240),
                  pw.SizedBox(width: 10),
                  pw.Expanded(flex: 1, child: pepRightTable),
                ],
              ),
            ],
          );
        },
      ),
    );

    String raTitle = reportTitle.toLowerCase().trim();
    final file = File("$outputPath/$raTitle.pdf");
    await file.writeAsBytes(await document.save());
  }
}

/// A model representation of a RNA sequence.
class RNA extends Nucleotides {
  /// Creates a `RNA` object.
  RNA({required String seq}) : super._(seq: seq, type: 'rna');

  /// Returns the reverse-transcribed version of this sequence.
  String revTranscribe() => this.seq.replaceAll('U', 'T');

  /// Returns a `RNA` object with a specified length of [len].
  static RNA random({required int len}) {
    Random _rand = Random();
    String rnaNucsStr = rnaNucs.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => rnaNucsStr.codeUnitAt(
          _rand.nextInt(rnaNucsStr.length),
        ),
      ),
    );
    return RNA(seq: seq);
  }
}

/// A model representation of a peptide sequence.
class Peptide extends Sequence {
  Peptide({required String seq}) : super._(seq: seq, type: 'pep');

  /// Returns the monoisotopic mass (Da) of this sequence.
  ///
  /// Use [roundTo] to specify the number of decimal places.
  /// Return the result in units of kDa by setting [kDa] to `true`.
  double monoMass({int roundTo = 3, kDa = false}) {
    double totalMonoMass = 0;
    for (var aa in this.seq.split('')) {
      totalMonoMass += aaToMonoMass[aa]!;
    }
    if (kDa) {
      return double.parse((totalMonoMass / 1000).toStringAsFixed(roundTo));
    }
    return double.parse(totalMonoMass.toStringAsFixed(roundTo));
  }

  /// Returns a `Peptide` object with a specified length of [len].
  static Peptide random({required int len}) {
    Random _rand = Random();
    String aminoAcidsStr = aminoAcids.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => aminoAcidsStr.codeUnitAt(
          _rand.nextInt(aminoAcidsStr.length),
        ),
      ),
    );
    return Peptide(seq: seq);
  }
}

// void main() async {
//   final fastaData = await Utils.readFASTA(path: '../example_dna_fasta.txt');
//   DNA dna = new DNA(seq: fastaData.first['seq']!);
//   dna.toFASTA(outputPath: '../deliverables', filename: 'fasta_output');
//   dna.report(outputPath: '../', creatorName: 'John Doe', reportTitle: 'Example Report');
// }
