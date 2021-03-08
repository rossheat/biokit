library biokit;

import 'dart:collection';

import 'package:pdf/pdf.dart';
import 'package:pdf/widgets.dart' as pw;
import 'dart:convert';
import 'dart:io';
import 'dart:math';

const String kDNA = 'dna';
const String kRNA = 'rna';
const String kPep = 'pep';
const String kAASeq = 'aaSeq';
const String kStartIndex = 'startIndex';
const String kEndIndex = 'endIndex';
const String kMatchCount = 'matchCount';
const String kMatchIndices = 'matchIndices';
const String kMatch = 'match';
const String kSeq = 'seq';
const String kId = 'id';
const String kDesc = 'desc';

class Structs {
  static const List<String> dnaNucs = ["A", "T", "G", "C"];
  static const List<String> rnaNucs = ["A", "U", "G", "C"];
  static const List<String> translationStopCodons = ["UGA", "UAA", "UAG"];

  static const List<String> aminoAcids = [
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

  static const List<String> dnaTransitions = ['CT', 'TC', 'AG', 'GA'];
  static const List<String> dnaTransversions = ['GT', 'TG', 'AC', 'CA', 'AT', 'TA', 'GC', 'CG'];

  static const Map<String, String> dnaCompNucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'};
  static const Map<String, String> rnaCompNucs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'};

  // X = STOP codon
  static const Map<String, String> dnaCodonToAA = {
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

  // X = STOP codon
  static const Map<String, String> rnaCodonToAA = {
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

  // Monoisotopic Mass
  static const Map<String, double> aaToMonoMass = {
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
}

class Errors {
  static String invalidSeq({required String mon, required int idx, required String type}) {
    return "Invalid ${type.toUpperCase()} Sequence Error. Character '$mon' found at index position $idx (zero-based) is not a valid ${type.toUpperCase()} monomer.";
  }

  static String invalidType({required String type}) {
    return "Invalid Sequence Type Error. '$type' is not a valid sequence type. Please enter the argument 'dna', 'rna' or 'pep' for [type].";
  }
}

class Utils {
  static Future<String> uniprotIdToFASTA({required String uniprotId}) async {
    Uri uri = Uri.parse('http://www.uniprot.org/uniprot/$uniprotId.fasta');

    var request = await HttpClient().getUrl(uri);
    var response = await request.close();

    await for (var contents in response.transform(Utf8Decoder())) {
      return contents;
    }
    return 'Error retrieving protein with uniprot ID $uniprotId';
  }

  static Future<List<Map<String, String>>> readFASTA({required String path}) async {
    String contents = await File(path).readAsString();
    List<String> lines = contents.split('\n');
    int seqCount = 0;

    List<Map<String, String>> fastaMaps = [];
    Map<String, String> currentMap = {};

    for (var line in lines) {
      if (line.startsWith('>')) {
        // Starting new line
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

  static String motifToRe({required motif}) {
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

class Sequence {
  late final String _seq;
  late final String _type;
  late final int _len;
  late String _name;
  late String _id;
  late String _desc;

  Sequence._({required String seq, required String type}) {
    this._type = _validateType(type: type);
    this._seq = _validateSeq(seq: seq);
    this._name = 'Default name';
    this._id = 'Default ID';
    this._desc = 'Default description';
  }

  String _validateType({required String type}) {
    String lType = type.toLowerCase();
    if (![kDNA, kRNA, kPep].contains(lType)) {
      throw (Errors.invalidType(type: lType));
    }
    return lType;
  }

  String _validateSeq({required String seq}) {
    int seqLen = seq.length;
    if (seqLen < 6) {
      throw ('Invalid Sequence Length Error. Sequence must have 6 or more monomers.');
    }

    String uSeq = seq.toUpperCase();
    uSeq.split('').asMap().forEach((idx, mon) {
      if (this._type == kDNA) {
        if (!Structs.dnaNucs.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      } else if (this.type == kRNA) {
        if (!Structs.rnaNucs.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      } else {
        if (!Structs.aminoAcids.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      }
    });
    this._len = seqLen;
    return uSeq;
  }

  String operator +(Sequence oSeq) {
    if (this._type != oSeq.type) {
      throw ('Cannot add ${this._type} and ${oSeq.type} sequence.');
    }
    return this._seq + oSeq.seq;
  }

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

  @override
  String toString() => info().toString();

  String get seq => this._seq;

  String get type => this._type;

  int get len => this._len;

  String get name => this._name;

  String get id => this._id;

  String get desc => this._desc;

  set name(String newName) {
    if (newName.length < 2 || 25 > newName.length) {
      throw ('Name Length Error. Name must between 2 and 25 characters.');
    }
    this._name = newName;
  }

  set id(String newId) {
    if (newId.length < 2 || 30 > newId.length) {
      throw ('ID Length Error. ID must between 2 and 30 characters.');
    }
    this._id = newId;
  }

  set desc(String newDesc) {
    if (newDesc.length < 5 || 100 > newDesc.length) {
      throw ('Description Length Error. Description must between 5 and 100 characters.');
    }
    this._desc = newDesc;
  }

  /// Calculates the frequency of each monomer.
  Map<String, double> freq({bool norm = false, bool ignoreStopCodon = true}) {
    Map<String, double> freqMap = {};
    this._seq.split('').forEach((mon) {
      if (this is Peptide && ignoreStopCodon) {
        if (mon != 'X') {
          if (freqMap.containsKey(mon)) {
            freqMap[mon] = freqMap[mon]! + 1;
          } else {
            freqMap[mon] = 1;
          }
        }
      } else {
        if (freqMap.containsKey(mon)) {
          freqMap[mon] = freqMap[mon]! + 1;
        } else {
          freqMap[mon] = 1;
        }
      }
    });

    if (ignoreStopCodon & (this is Peptide)) {
      if (norm) {
        return freqMap.map(
          (key, value) => MapEntry(
            key,
            double.parse(((value / _lenMinusStopAA(stopAA: 'X')) * 100).toStringAsFixed(1)),
          ),
        );
      }
      return freqMap;
    } else {
      if (norm) {
        return freqMap.map(
          (key, value) => MapEntry(
            key,
            double.parse(((value / this._len) * 100).toStringAsFixed(1)),
          ),
        );
      }
      return freqMap;
    }
  }

  String _reversed({required String seq}) => seq.split('').reversed.join('');

  /// Reverses the sequence.
  String reverse() => this._seq.split('').reversed.join('');

  /// Every possible comnbination of the sequence.
  List<String> combinations({sorted = false}) {
    List<String> listSeq = this._seq.split("");
    List<String> combinations = [];
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
      // Sort with longest combination first.
      combinations.sort((b, a) => a.length.compareTo(b.length));
      return combinations;
    }
    return combinations;
  }

  /// Finds the indices (zero-based) of a specified motif.
  Map<String, dynamic> findMotif({required String motif, overlap = true}) {
    List<Map<String, dynamic>> matchData = [];
    Map<String, dynamic> matchMotifMap = {};
    String tempRegexMotif = Utils.motifToRe(motif: motif);
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

  // The number of positional differences.
  int difference({required Sequence oSeq}) {
    if (this._len != oSeq.len) {
      throw ('Sequences must be of the same length to calculate difference.');
    }
    if (this._type != oSeq.type) {
      throw ('Sequences must be of the same type to calculate difference.');
    }

    int differenceCount = 0;
    this._seq.split('').asMap().forEach((idx, mon) {
      if (mon != oSeq.seq[idx]) {
        differenceCount++;
      }
    });
    return differenceCount;
  }

  /// Removes all occurrences of the specified motif.
  String splice({required String motif}) {
    String vMotif = _validateSeq(seq: motif);
    return seq.replaceAll(vMotif, '');
  }

  // The longest shared motif.
  String sharedMotif({required oSeq}) {
    if (this._type != oSeq.type) {
      throw ('Cannot find shared motif between ${this._type} and ${oSeq.type} sequence.');
    }

    // Generate all possible combinations.
    List<String> combos = combinations(sorted: true);

    String longestShared = '';

    // Find the longest combination that is contained in all sequences.
    for (var comb in combos) {
      bool allMatches = true;
      if (!oSeq.seq.contains(comb)) {
        allMatches = false;
        break;
      }
      if (allMatches) {
        longestShared = comb;
        break;
      }
    }
    return longestShared;
  }

  double _norm({required int number, required int total, required int places}) =>
      double.parse(((number / total) * 100).toStringAsFixed(1));

  int _lenMinusStopAA({required String stopAA}) {
    int tempLen = 0;
    for (var aa in this._seq.split('')) {
      if (aa != stopAA) {
        tempLen++;
      }
    }
    return tempLen;
  }
}

class Nucleotides extends Sequence {
  Nucleotides._({required String seq, required String type}) : super._(seq: seq, type: type);

  // /// The frequency of each nucleotide.
  // Map<String, double> freq({bool norm = false}) => super.freq(norm: norm);

  /// Translates the nucleotides to amino acids.
  Map<String, dynamic> translate({revComp = false, startIdx = 0}) {
    String seq = revComp ? complementary(rev: true) : this.seq;

    String aaSeq = '';
    for (var i = startIdx; i < seq.length - 2; i += 3) {
      String codon = seq.substring(i, i + 3);
      aaSeq += this.type == kDNA ? Structs.dnaCodonToAA[codon]! : Structs.rnaCodonToAA[codon]!;
    }
    return {kAASeq: aaSeq, 'nucCount': seq.length - startIdx - 1, 'aaCount': aaSeq.length};
  }

  /// Calculates the frequency of each codon for a specified amino acid.
  /// Searches in batches of three nucleotides.
  /// If a codon is present but does not appear in a batch, it will not be detected.
  /// For example in ATGTCATGC, only 1 ATG is detected.
  // Map<String, int> codonFreq({required String aminoAcid}) {
  //   Map<String, int> codonFreqMap = {};
  //   for (var i = 0; i < this._len - 2; i += 3) {
  //     String codon = this.seq.substring(i, i + 3);
  //     String fetchedAminoAcid =
  //         this.type == kDNA ? Structs.dnaCodonToAA[codon]! : Structs.rnaCodonToAA[codon]!;
  //     if (fetchedAminoAcid == aminoAcid.toUpperCase()) {
  //       codonFreqMap[codon] == null
  //           ? codonFreqMap[codon] = 1
  //           : codonFreqMap[codon] = codonFreqMap[codon]! + 1;
  //     }
  //   }
  //   return codonFreqMap;
  // }

  int codonFreq({required String codon}) {
    int codonFreq = 0;
    for (var i = 0; i < this._len - 2; i += 3) {
      String fetchedCodon = this._seq.substring(i, i + 3);
      if (codon == fetchedCodon) {
        codonFreq++;
      }
    }
    return codonFreq;
  }

  /// The complementary strand.
  String complementary({bool rev = false}) {
    String compSeq = this
        .seq
        .split('')
        .map((nuc) => this.type == kDNA ? Structs.dnaCompNucs[nuc] : Structs.rnaCompNucs[nuc])
        .join();
    return rev ? super._reversed(seq: compSeq) : compSeq;
  }

  /// The percentage of nucleotides which are either Guanine or Cytosine.
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

  /// Generates 6 ORFs, 3 forward, 3 reversed.
  List<String> readingFrames() {
    List<String> readingFrames = [];

    for (var i = 0; i < 3; i++) {
      readingFrames.add(translate(startIdx: i)[kAASeq]);
      readingFrames.add(translate(revComp: true, startIdx: i)[kAASeq]);
    }
    return readingFrames;
  }

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

  List<String> proteins({bool unique = false}) {
    List<String> frames = readingFrames();
    List<String> proteins = [];
    for (var frame in frames) {
      List<String> tempProteins = readingFrameToProteins(aaSeq: frame);
      if (proteins != []) {
        proteins.addAll(tempProteins);
      }
    }
    proteins.sort((b, a) => a.length.compareTo(b.length));
    if (unique) {
      return proteins.toSet().toList();
    }
    return proteins;
  }

  double molWeight() => this._len * 0.33; // kDa

  double doubleHelixMolWeight() => molWeight() * 2; // kDa

  double doubleHexlixTurns() => this._len / 10; // turns

  double doubleHelixGeoLen() => this._len * 0.34; //nm
}

class DNA extends Nucleotides {
  DNA({required String seq}) : super._(seq: seq, type: 'dna');

  /// Transcribes DNA to RNA.
  String transcribe() => this.seq.replaceAll('T', 'U');

  Map<String, List<Map<String, int>>> restrictionSites({
    int minSiteLen = 4,
    int maxSiteLen = 8,
  }) {
    String revComp = complementary(rev: true);

    List<String> origSeqCombos = combinations();
    Iterable<String> seqIter =
        origSeqCombos.where((seq) => (seq.length >= minSiteLen) && (seq.length <= maxSiteLen));
    List<String> restSeqs = seqIter.toSet().toList();

    Map<String, List<Map<String, int>>> restSiteSeqs = {};
    List<Map<String, int>> restSiteLocations = [];

    for (String restSeq in restSeqs) {
      restSiteLocations = [];
      RegExp regRestExp = RegExp(restSeq);
      Iterable<RegExpMatch> matches = regRestExp.allMatches(revComp);
      for (RegExpMatch match in matches) {
        int startIdx = this.seq.length - (restSeq.length + match.start);
        int endIdx = startIdx + restSeq.length;
        if (this.seq.substring(startIdx, endIdx) == restSeq) {
          restSiteLocations.add({kStartIndex: startIdx, kEndIndex: endIdx});
          restSiteSeqs[restSeq] = restSiteLocations;
        }
        ;
      }
    }
    return restSiteSeqs;
  }

  double tranRatio({required DNA oSeq}) {
    if (this._len != oSeq._len) {
      throw ('Unequal Sequence Lengths Error.');
    }

    int transitionCount = 0;
    int transversionCount = 0;

    this.seq.split('').asMap().forEach((idx, nuc) {
      if (nuc != oSeq.seq[idx]) {
        if (Structs.dnaTransitions.contains(nuc + oSeq.seq[idx])) {
          transitionCount++;
        } else if (Structs.dnaTransversions.contains(nuc + oSeq.seq[idx])) {
          transversionCount++;
        }
      }
    });
    return transitionCount / transversionCount;
  }

  static DNA random({required int len}) {
    Random _rand = Random();
    String dnaNucsStr = Structs.dnaNucs.join();
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

  Future<void> report(
      {required String outputPath,
      required String creatorName,
      required String reportTitle}) async {
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
      pep._lenMinusStopAA(stopAA: 'X'),
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

class RNA extends Nucleotides {
  RNA({required String seq}) : super._(seq: seq, type: 'rna');

  /// Transcribe RNA back to DNA.
  String revTranscribe() => this.seq.replaceAll('U', 'T');

  static RNA random({required int len}) {
    Random _rand = Random();
    String rnaNucsStr = Structs.rnaNucs.join();
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

class Peptide extends Sequence {
  Peptide({required String seq}) : super._(seq: seq, type: 'pep');

  /// Calculate the Monoisotopic Mass.
  double monoMass({int roundTo = 3, kDa = false}) {
    double totalMonoMass = 0;
    for (var aa in this.seq.split('')) {
      totalMonoMass += Structs.aaToMonoMass[aa]!;
    }
    if (kDa) {
      return double.parse((totalMonoMass / 1000).toStringAsFixed(roundTo));
    }
    return double.parse(totalMonoMass.toStringAsFixed(roundTo));
  }

  static Peptide random({required int len}) {
    Random _rand = Random();
    String aminoAcidsStr = Structs.aminoAcids.join();
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

void main() async {
  // List<Map<String, String>> dnaSeqs = await Utils.readFASTA(path: 'example_dna_fasta.txt');
  // DNA dna = DNA(seq: dnaSeqs.first['seq']!);
  // dna.report(
  //   outputPath: '../deliverables',
  //   reportTitle: 'DNA Analysis Report',
  //   creatorName: 'John Doe',
  // );
}
