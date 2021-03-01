import 'package:biokit/biokit.dart';
import 'package:flutter/material.dart';

void main() {
  runApp(MyApp());
}

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'BioKit Demo',
      home: HomeScreen(),
    );
  }
}

class HomeScreen extends StatefulWidget {
  @override
  _HomeScreenState createState() => _HomeScreenState();
}

class _HomeScreenState extends State<HomeScreen> {
  @override
  void initState() {
    super.initState();

    String DNASequence = makeNucleotideSequence(40, 'DNA');
    String RNASequence = makeNucleotideSequence(40, 'RNA');
    print(validateNucleotideSequence(DNASequence, 'DNA'));
    print(nucleotideFrequency(DNASequence, 'DNA'));
    print(DNAToRNA('GATGGAACTTGACTACGTAAATT'));
    print(complementaryStrand('AAAACCCGGT', 'DNA', true));
    print(GCContent('CCACCCTCGTGGTAGGCAGTAGGTGGAAT', 'DNA'));
    print(hammingDistance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT', 'DNA'));
    print(dominantPhenotypeProba(19, 29, 28)); // 0.6892982456140351
    print(detectSequenceType('AGCGG'));
    print(RNAtoProtein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'));
    print(matchMotif('GATATATGCATATACTT', 'ATAT', 'DNA'));
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      body: Center(
        child: Text('BioKit'),
      ),
    );
  }
}
