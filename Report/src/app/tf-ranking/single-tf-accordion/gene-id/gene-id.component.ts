import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'data-gene-id',
  templateUrl: './gene-id.component.html',
  styleUrls: ['./gene-id.component.css']
})
export class GeneIDComponent implements OnInit {
  @Input()
    // @ts-ignore
  geneID: string;

  constructor() {
  }

  ngOnInit(): void {
  }

}
