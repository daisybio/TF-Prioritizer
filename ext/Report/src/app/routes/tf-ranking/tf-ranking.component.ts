import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";
import {TfDataGetterService} from "../../services/tf-data-getter.service";
import {resourceChangeTicket} from "@angular/compiler-cli/src/ngtsc/core";

@Component({
  selector: 'tf-ranking',
  templateUrl: './tf-ranking.component.html',
  styleUrls: ['./tf-ranking.component.css']
})
export class TfRankingComponent implements OnInit {
  @Input()
  tfGroups: TranscriptionFactorGroup[];

  constructor(tfDataGetter: TfDataGetterService) {
    this.tfGroups = tfDataGetter.getData().transcriptionFactorGroups;
  }

  ngOnInit(): void {
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
