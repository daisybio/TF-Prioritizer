import {Injectable} from '@angular/core';
import {HttpClient} from "@angular/common/http";
import {ranking} from "../../interfaces";
import * as ranks from 'src/assets/data/ranks.json';
import {TranscriptionFactor} from "../classes/transcription-factor";

@Injectable({
  providedIn: 'root'
})
export class DataService {
  ranking: ranking;
  transcriptionFactors: Set<TranscriptionFactor> = new Set<TranscriptionFactor>();

  constructor(private httpClient: HttpClient) {
    this.ranking = ranks;

    let allTfNames = new Set<string>();

    for (let assay of this.getAssays()) {
      for (let pairing of this.getPairings()) {

        if (!(pairing in this.ranking[assay])) {
          continue;
        }

        for (let tf of this.ranking[assay][pairing]) {
          allTfNames.add(tf);
        }
      }
    }

    for (let tfName of allTfNames) {
      this.transcriptionFactors.add(new TranscriptionFactor(tfName, this.ranking));
    }
  }

  public getRanking(): ranking {
    return this.ranking
  }

  public getAssays(): string[] {
    return Object.keys(this.ranking).filter(val => val != 'default').sort();
  }

  public getPairings(): string[] {
    let assays = this.getAssays();

    let pairings = new Set<string>();

    for (let assay of assays) {
      for (let pairing of Object.keys(this.ranking[assay])) {
        pairings.add(pairing);
      }
    }

    return Array.from(pairings).sort();
  }
}
