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
    String DNASequence = generateRandomNucleotideSequence(40, 'DNA');
    String RNASequence = generateRandomNucleotideSequence(40, 'RNA');
    print(validateNucleotideSequence(DNASequence, 'DNA'));
    print(countNucleotideFrequency(DNASequence, 'DNA'));
    print(transcribeDNAToRNA('GATGGAACTTGACTACGTAAATT'));
    print(generateComplementaryStrand('AAAACCCGGT', 'DNA', true));
    print(calculateGCContent('CCACCCTCGTGGTAGGCAGTAGGTGGAAT', 'DNA'));
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
