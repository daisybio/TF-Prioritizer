import {Component, EventEmitter, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'tf-ranking',
  templateUrl: './tf-ranking.component.html',
  styleUrls: ['./tf-ranking.component.css']
})
export class TfRankingComponent implements OnInit {
  tfGroup_rank: { [tfGroup: string]: number } = {};
  existingHms = this.tfDataGetter.getData().existingValues.hm;
  activeHms: Set<string> = new Set<string>();
  rankingChangeEmitter: EventEmitter<{ [tfGroup: string]: number }> = new EventEmitter();

  constructor(private tfDataGetter: TfDataGetterService) {
  }

  ngOnInit(): void {
    this.existingHms.forEach(hm => this.activeHms.add(hm));
    this.updateRanking();
  }

  getTfGroups() {
    return Object.keys(this.tfGroup_rank);
  }

  updateRanking() {
    let ranked = this.tfDataGetter.getRanked(Array.from(this.activeHms.values()));
    this.tfGroup_rank = {};
    for (let i = 0; i < ranked.length; i++) {
      this.tfGroup_rank[ranked[i]] = i;
    }
    this.rankingChangeEmitter.emit(this.tfGroup_rank);
  }

  isRanked(tfGroup: string) {
    return Object.keys(this.tfGroup_rank).includes(tfGroup);
  }

  getRank(tfGroup: string) {
    return this.tfGroup_rank[tfGroup];
  }

  getTfGroup(name: string) {
    return this.tfDataGetter.getTfGroupByName(name);
  }

  toggleHm(hm: string) {
    if (this.activeHms.has(hm)) {
      if (!this.isOnlyActiveHm(hm))
        this.activeHms.delete(hm);
    } else {
      this.activeHms.add(hm);
    }
    this.updateRanking();
  }

  isOnlyActiveHm(hm: string) {
    return this.activeHms.size === 1 && this.activeHms.has(hm);
  }

  hasTargetGene(term: string, tfGroup: TranscriptionFactorGroup) {
    for (let targetGene of tfGroup.targetGenes) {
      if (targetGene.geneID.toUpperCase().includes(term.toUpperCase()) || targetGene.symbol.toUpperCase().includes(term.toUpperCase())) {
        return true;
      }
    }
    return false;
  }
}
