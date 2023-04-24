import {Component, OnInit} from '@angular/core';
import {DataService} from "../../services/data.service";
import {Subject} from "rxjs";
import {TranscriptionFactor} from "../../classes/transcription-factor";

@Component({
  selector: 'app-tf-ranking',
  templateUrl: './tf-ranking.component.html',
  styleUrls: ['./tf-ranking.component.scss']
})
export class TfRankingComponent implements OnInit {
  exisingAssays: string[];
  exisingPairings: string[];
  transcriptionFactors: Set<TranscriptionFactor>;
  $ranking: Subject<TranscriptionFactor[]> = new Subject<TranscriptionFactor[]>();
  $activePairings = new Subject<Set<string>>;
  $activeAssays = new Subject<Set<string>>;

  activePairings: Set<string> = new Set<string>();
  activeAssays: Set<string> = new Set<string>();

  constructor(private data: DataService) {
    this.exisingAssays = this.data.getAssays();
    this.exisingPairings = this.data.getPairings();
    this.transcriptionFactors = data.transcriptionFactors;

    this.$activeAssays.subscribe(val => {
      this.activeAssays = val;
      this.updateRanking();
    });

    this.$activePairings.subscribe(val => {
      this.activePairings = val;
      this.updateRanking();
    });

    this.$activePairings.next(new Set<string>(this.exisingPairings));
    this.$activeAssays.next(new Set<string>(this.exisingAssays));
  }

  ngOnInit(): void {
  }

  ngAfterViewInit(): void {
    setTimeout(() => {
      this.updateRanking();
    });
  }

  private updateRanking() {
    let ranking = Array.from(this.transcriptionFactors).map(tf => {
      return {
        tf: tf,
        dcg: tf.getDcgSum(this.activeAssays, this.activePairings)
      }
    }).filter(val => val.dcg > 0)
      .sort((a, b) => b.dcg - a.dcg)
      .map(val => val.tf);

    this.$ranking.next(ranking);
  }
}
