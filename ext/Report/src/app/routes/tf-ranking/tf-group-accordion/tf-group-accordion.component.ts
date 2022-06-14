import {Component, EventEmitter, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../../types/types";

@Component({
  selector: 'tf-group-accordion',
  templateUrl: './tf-group-accordion.component.html',
  styleUrls: ['./tf-group-accordion.component.css']
})
export class TfGroupAccordionComponent implements OnInit {
  @Input()
  tfGroup!: TranscriptionFactorGroup;

  @Input()
  showAnalysisLinks: boolean = true;

  @Input()
  index: number | undefined;
  @Input()
  rankingChangeEmitter?: EventEmitter<{ [tfGroup: string]: number }>

  visible: boolean = false
  title: string = "Undefined";

  constructor() {
  }

  ngOnInit(): void {
    this.updateTitle();

    this.rankingChangeEmitter?.subscribe((map: { [tfGroup: string]: number }) => {
      this.index = map[this.tfGroup.name];
      this.updateTitle();
    });
  }

  updateTitle() {
    if (this.index === undefined) {
      this.title = "Basic data";
    } else {
      this.title = (this.index + 1) + ". " + this.tfGroup.name;
    }
  }
}
